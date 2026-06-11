"""Differential expression and abundance across sample groups.

single_cell: Wilcoxon on pooled cells (high power, pseudoreplication).
pseudobulk:  per-sample means compared across samples (correct between-sample
             test, conservative with few samples).
Fold-changes are computed from linear raw intensities.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

from .config import PipelineConfig


def _bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _log2fc_linear(adata, markers, mask_a, mask_b, eps=1.0):
    raw = adata[:, markers].layers["raw"]
    ma = raw[mask_a].mean(axis=0)
    mb = raw[mask_b].mean(axis=0)
    return np.log2((ma + eps) / (mb + eps)), ma, mb


def single_cell_dea(adata, groupby, group, cfg: PipelineConfig):
    """One-vs-rest Wilcoxon on pooled cells for a single group value."""
    markers = adata.uns["channels"]["markers"]
    sub = adata[:, markers].copy()
    sub.X = sub.layers["zscore"]
    sub.obs["_grp"] = np.where(sub.obs[groupby].astype(str) == group, group, "rest")
    sub.obs["_grp"] = sub.obs["_grp"].astype("category")
    if sub.obs["_grp"].value_counts().min() < cfg.min_cells_per_group:
        return None
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sc.tl.rank_genes_groups(sub, "_grp", groups=[group], reference="rest",
                                method=cfg.dea_method)
    res = sc.get.rank_genes_groups_df(sub, group=group)
    mask_a = (sub.obs["_grp"] == group).to_numpy()
    log2fc, ma, mb = _log2fc_linear(adata, markers, mask_a, ~mask_a)
    fc = pd.Series(log2fc, index=markers)
    res = res.rename(columns={"names": "marker", "pvals": "pval",
                              "pvals_adj": "padj", "scores": "wilcoxon_z"})
    res["log2FC"] = res["marker"].map(fc)
    res["mean_group"] = res["marker"].map(pd.Series(ma, index=markers))
    res["mean_rest"] = res["marker"].map(pd.Series(mb, index=markers))
    res["group"] = group
    res["neglog10_padj"] = -np.log10(res["padj"].clip(lower=1e-300))
    return res.sort_values("log2FC", ascending=False).reset_index(drop=True)


def pseudobulk_dea(adata, groupby, sample_key="sample"):
    """Per-sample mean per marker, then Mann-Whitney U across samples.

    Only defined for a two-group comparison. Returns tidy df comparing the two
    largest groups (group_a vs group_b) found in `groupby`.
    """
    markers = adata.uns["channels"]["markers"]
    df = pd.DataFrame(adata[:, markers].layers["transformed"], columns=markers,
                      index=adata.obs_names)
    df[sample_key] = adata.obs[sample_key].values
    df["_grp"] = adata.obs[groupby].astype(str).values
    # map sample -> its (single) group
    samp_grp = df.groupby(sample_key, observed=True)["_grp"].agg(
        lambda s: s.value_counts().index[0])
    pb = df.groupby(sample_key, observed=True)[markers].mean()
    groups = samp_grp.value_counts().index.tolist()
    if len(groups) < 2:
        return None, None
    # use a healthy/control group as the reference where one is present
    ref_tokens = ("healthy", "control", "normal", "donor", "hd", "ctrl", "wt")
    ctrl = [g for g in groups if any(t in g.lower() for t in ref_tokens)]
    if ctrl:
        gb = ctrl[0]
        ga = next(g for g in groups if g != gb)
    else:
        ga, gb = sorted(groups)[0], sorted(groups)[1]
    a_samples = samp_grp[samp_grp == ga].index
    b_samples = samp_grp[samp_grp == gb].index

    rawdf = pd.DataFrame(adata[:, markers].layers["raw"], columns=markers,
                         index=adata.obs_names)
    rawdf[sample_key] = adata.obs[sample_key].values
    pb_raw = rawdf.groupby(sample_key, observed=True)[markers].mean()

    rows = []
    for m in markers:
        a = pb.loc[a_samples, m].to_numpy()
        b = pb.loc[b_samples, m].to_numpy()
        try:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            p = 1.0
        ra = pb_raw.loc[a_samples, m].mean()
        rb = pb_raw.loc[b_samples, m].mean()
        log2fc = np.log2((ra + 1.0) / (rb + 1.0))
        # standardized effect (Hedges-ish) on transformed means
        pooled = np.sqrt(((a.var(ddof=1) if len(a) > 1 else 0) +
                          (b.var(ddof=1) if len(b) > 1 else 0)) / 2 + 1e-9)
        d = (a.mean() - b.mean()) / pooled if pooled > 0 else 0.0
        rows.append((m, log2fc, d, p, ra, rb))
    out = pd.DataFrame(rows, columns=["marker", "log2FC", "cohens_d", "pval",
                                      "mean_group", "mean_rest"])
    out["padj"] = _bh(out["pval"].to_numpy())
    out["neglog10_padj"] = -np.log10(out["padj"].clip(lower=1e-300))
    out["group"] = ga
    out["reference"] = gb
    meta = {"group_a": ga, "group_b": gb,
            "n_a": len(a_samples), "n_b": len(b_samples),
            "samples_a": list(a_samples), "samples_b": list(b_samples)}
    return out.sort_values("log2FC", ascending=False).reset_index(drop=True), meta


def differential_abundance(adata, groupby, cluster_key="cell_type_auto"):
    """Fraction of each group's cells in each cluster/cell type (col-normalized)."""
    ct = pd.crosstab(adata.obs[cluster_key], adata.obs[groupby])
    frac = ct.div(ct.sum(axis=0).replace(0, 1), axis=1)
    return ct, frac


def run_dea(adata, cfg: PipelineConfig):
    groupby = cfg.dea_groupby
    if groupby not in adata.obs or adata.obs[groupby].astype(str).nunique() < 2:
        return None
    groups = [g for g in adata.obs[groupby].astype(str).unique() if g != ""]
    result = {"groupby": groupby, "groups": groups, "single_cell": {}}
    for g in groups:
        r = single_cell_dea(adata, groupby, g, cfg)
        if r is not None:
            result["single_cell"][g] = r
    pb, meta = pseudobulk_dea(adata, groupby)
    result["pseudobulk"] = pb
    result["pseudobulk_meta"] = meta
    abund_key = "cell_type_auto" if "cell_type_auto" in adata.obs else "leiden"
    ct, frac = differential_abundance(adata, groupby, abund_key)
    result["abundance_counts"] = ct
    result["abundance_frac"] = frac
    result["abundance_key"] = abund_key
    return result
