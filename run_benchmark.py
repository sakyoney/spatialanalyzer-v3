"""Benchmark the MACSima Analyzer software against synthetic ground truth.

This is a SEPARATE step from the software. It:
  1. generates a synthetic Multiple Myeloma dataset with known ground truth,
  2. runs the (unmodified) software pipeline on it,
  3. scores how well the pipeline recovered the seeded biology,
  4. writes an example interactive report.

The software package (`macsima`) has no knowledge of any of this - the
dependency points one way only: benchmark -> software.

Run from the project root:
    python benchmark/run_benchmark.py
"""
from __future__ import annotations
import os
import sys
import time

# make the software package importable, and this folder for generate_mm
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)
sys.path.insert(0, _HERE)

import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from macsima.config import PipelineConfig
from macsima import pipeline, report
import generate_mm

# Markers the generator deliberately UP-regulates on myeloma plasma cells.
SEEDED_DISEASE_MARKERS = ["BCMA", "GPRC5D", "CD138", "CD38"]


def _match(markers, token):
    return next((m for m in markers if token.lower() in m.lower()), None)


def score_annotation(adata):
    true = adata.obs["_true_type"].astype(str).to_numpy()
    auto = adata.obs["cell_type_auto"].astype(str).to_numpy()
    ari = adjusted_rand_score(true, auto)
    nmi = normalized_mutual_info_score(true, auto)
    # plasma-cell recovery (the clinically important population)
    true_pc = (np.char.find(np.char.lower(true.astype(str)), "plasma") >= 0)
    auto_pc = (np.char.find(np.char.lower(auto.astype(str)), "plasma") >= 0)
    tp = int((true_pc & auto_pc).sum())
    precision = tp / max(auto_pc.sum(), 1)
    recall = tp / max(true_pc.sum(), 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-9)
    return {"ARI": ari, "NMI": nmi, "n_true_plasma": int(true_pc.sum()),
            "n_auto_plasma": int(auto_pc.sum()), "plasma_precision": precision,
            "plasma_recall": recall, "plasma_F1": f1}


def score_dea(results, top_k=6):
    dea = results["dea"]
    disease = next((g for g in dea["single_cell"] if "myeloma" in g.lower()),
                   list(dea["single_cell"])[0])
    df = dea["single_cell"][disease]
    top = df.nlargest(top_k, "log2FC")["marker"].tolist()
    hits = [tok for tok in SEEDED_DISEASE_MARKERS
            if any(tok.lower() in m.lower() for m in top)]
    return {"group": disease, "top_markers": top,
            "seeded_recovered": hits, "recovery": len(hits) / len(SEEDED_DISEASE_MARKERS)}


def score_spatial(results, top_k=6):
    moran = results["spatial"].get("moran")
    if moran is None or not len(moran):
        return {"recovery": float("nan"), "top_markers": []}
    top = moran.sort_values("morans_I", ascending=False).head(top_k)["marker"].tolist()
    hits = [tok for tok in SEEDED_DISEASE_MARKERS
            if any(tok.lower() in m.lower() for m in top)]
    return {"top_markers": top, "seeded_recovered": hits,
            "recovery": len(hits) / len(SEEDED_DISEASE_MARKERS)}


def main():
    data_dir = os.path.join(_HERE, "data", "synthetic_mm")
    print("=" * 70)
    print("MACSima Analyzer - benchmark")
    print("=" * 70)
    print(f"\n[1/4] Generating synthetic dataset -> {data_dir}")
    generate_mm.generate_dataset(data_dir)

    print("\n[2/4] Running software pipeline (macsima)...")
    cfg = PipelineConfig()
    t0 = time.time()
    adata, results = pipeline.run_pipeline(
        data_dir, cfg,
        progress=lambda s, m, dt: print(f"      [{s:9s}] {m} ({dt:.1f}s)"))
    runtime = time.time() - t0

    print("\n[3/4] Scoring recovery vs ground truth...")
    ann = score_annotation(adata)
    dea = score_dea(results)
    spat = score_spatial(results)

    print("\n" + "-" * 70)
    print("ANNOTATION (unsupervised auto vs true cell type)")
    print(f"  Adjusted Rand Index : {ann['ARI']:.3f}   (1.0 = perfect)")
    print(f"  Normalised MI       : {ann['NMI']:.3f}")
    print(f"  Plasma cells        : true={ann['n_true_plasma']}, "
          f"auto={ann['n_auto_plasma']}, F1={ann['plasma_F1']:.3f} "
          f"(P={ann['plasma_precision']:.3f}, R={ann['plasma_recall']:.3f})")
    print("\nDIFFERENTIAL EXPRESSION (top up-regulated in disease, single-cell)")
    print(f"  Comparison group    : {dea['group']}")
    print(f"  Top markers         : {', '.join(dea['top_markers'])}")
    print(f"  Seeded recovered    : {', '.join(dea['seeded_recovered'])} "
          f"({dea['recovery']*100:.0f}% of {len(SEEDED_DISEASE_MARKERS)})")
    print("\nSPATIAL AUTOCORRELATION (top Moran's I)")
    print(f"  Top markers         : {', '.join(spat['top_markers'])}")
    print(f"  Seeded recovered    : {', '.join(spat.get('seeded_recovered', []))} "
          f"({spat['recovery']*100:.0f}%)")
    print("-" * 70)

    overall = np.mean([ann["plasma_F1"], dea["recovery"], spat["recovery"]])
    verdict = "PASS" if overall >= 0.8 else "REVIEW"
    print(f"\nOverall recovery score: {overall:.2f}  ->  {verdict}")
    print(f"Pipeline runtime: {runtime:.1f}s  ({adata.n_obs:,} cells)")

    print("\n[4/4] Writing example report...")
    out = os.path.join(_HERE, "example_report.html")
    cfg.report_title = "MACSima Analyzer - Synthetic Benchmark Report"
    report.build_report(adata, results, cfg, out)
    print(f"      {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
    print("\nDone.")


if __name__ == "__main__":
    main()
