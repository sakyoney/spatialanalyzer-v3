"""MACSima Analyzer - Streamlit GUI.

Run with:   streamlit run app.py

Workflow: load data -> review channel classification -> tune each stage's
parameters (recommended defaults pre-filled) -> run -> download interactive
HTML report. Every box has a sensible default; you only touch what your
science requires.
"""
from __future__ import annotations
import os
import time
import tempfile
import streamlit as st

from macsima.config import PipelineConfig
from macsima import io, pipeline, report

st.set_page_config(page_title="MACSima Analyzer", layout="wide",
                   initial_sidebar_state="expanded")

STAGE_LABELS = dict(pipeline.STAGES)

# --------------------------------------------------------------------------
# session state
# --------------------------------------------------------------------------
ss = st.session_state
ss.setdefault("adata", None)
ss.setdefault("results", None)
ss.setdefault("stage_status", {})
ss.setdefault("loaded_adata", None)
ss.setdefault("csv_paths", [])
ss.setdefault("report_path", None)

st.title("MACSima Analyzer")
st.caption("Configurable spatial-proteomics pipeline for Miltenyi MACSima / "
           "MACS iQ View output - Scanpy + Squidpy + interactive HTML reports.")

# --------------------------------------------------------------------------
# Sidebar: data source
# --------------------------------------------------------------------------
with st.sidebar:
    st.header("1 | Data")
    source = st.radio("Source", ["Upload CSV files", "Folder path"])

    if source.startswith("Upload"):
        ups = st.file_uploader("MACSima CSVs", type=["csv"], accept_multiple_files=True)
        if ups:
            d = os.path.join(tempfile.gettempdir(), "macsima_upload")
            os.makedirs(d, exist_ok=True)
            paths = []
            for u in ups:
                p = os.path.join(d, u.name)
                with open(p, "wb") as f:
                    f.write(u.getbuffer())
                paths.append(p)
            ss.csv_paths = paths
            ss.loaded_adata = None
            st.success(f"{len(paths)} file(s) ready.")
    else:
        folder = st.text_input("Folder containing *.csv")
        if folder and os.path.isdir(folder):
            import glob
            ss.csv_paths = sorted(glob.glob(os.path.join(folder, "*.csv")))
            ss.loaded_adata = None
            st.success(f"{len(ss.csv_paths)} CSV(s) found.")

    st.divider()
    st.caption("Files queued:")
    for p in ss.csv_paths:
        st.write("-", os.path.basename(p))

cfg = PipelineConfig()

# --------------------------------------------------------------------------
# Load + channel classification preview
# --------------------------------------------------------------------------
if ss.csv_paths and ss.loaded_adata is None:
    try:
        ss.loaded_adata = io.load_folder(ss.csv_paths, cfg)
    except Exception as e:
        st.error(f"Import failed: {e}")

if ss.loaded_adata is not None:
    a = ss.loaded_adata
    st.subheader("2 | Channel classification")
    ch = a.uns["channels"]
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Cells", f"{a.n_obs:,}")
    m2.metric("Biological markers", len(ch["markers"]))
    m3.metric("DAPI / cycle", len(ch["dapi"]))
    m4.metric("Controls + blanks", len(ch["controls"]) + len(ch["bare_fluorophores"]))
    with st.expander("Review / override channel sets", expanded=False):
        st.write("**Markers (used downstream):**", ", ".join(ch["markers"]))
        st.write("**Excluded - DAPI:**", ", ".join(ch["dapi"]) or "-")
        st.write("**Excluded - bare fluorophores:**", ", ".join(ch["bare_fluorophores"]) or "-")
        st.write("**Excluded - controls/blanks:**", ", ".join(ch["controls"]) or "-")
        keep = st.multiselect("Force-include as markers", options=ch["dapi"] + ch["controls"] + ch["bare_fluorophores"])
        drop = st.multiselect("Force-exclude from markers", options=ch["markers"])

    # condition mapping
    st.subheader("3 | Sample -> condition mapping")
    samples = list(a.obs["sample"].cat.categories)
    existing = {s: a.obs.loc[a.obs["sample"] == s, "condition"].astype(str).iloc[0]
                for s in samples}
    cond_map = {}
    cols = st.columns(min(3, len(samples)))
    for i, s in enumerate(samples):
        default = existing.get(s, "")
        cond_map[s] = cols[i % len(cols)].text_input(f"{s}", value=default, key=f"cond_{s}")

# --------------------------------------------------------------------------
# Parameter boxes (workflow configuration)
# --------------------------------------------------------------------------
st.subheader("4 | Pipeline parameters")
st.caption("Defaults are recommended starting points - change only what your science requires.")
pc1, pc2, pc3 = st.columns(3)

with pc1:
    with st.expander("Transform & normalization", expanded=True):
        cfg.transform = st.selectbox("Transform", ["arcsinh", "log1p", "none"], 0)
        cfg.cofactor = st.number_input("arcsinh cofactor", 0.1, 1000.0, 5.0, 0.5)
        cfg.zscore = st.checkbox("Joint z-score (recommended)", True)
        st.caption("Joint = standardize after concatenation; preserves between-sample signal for DEA.")
    with st.expander("Clustering", expanded=True):
        cfg.use_pca = st.checkbox("Use PCA", True)
        cfg.n_pcs = st.slider("PCs", 2, 40, 15)
        cfg.n_neighbors = st.slider("k neighbors", 5, 50, 15)
        cfg.leiden_resolution = st.slider("Leiden resolution", 0.1, 3.0, 1.0, 0.1)

with pc2:
    with st.expander("Annotation", expanded=True):
        cfg.auto_annotate = st.checkbox("Auto cell-type annotation", True)
        cfg.positivity_z = st.slider("Positivity z-threshold", 0.0, 3.0, 1.0, 0.1)
    with st.expander("Differential expression", expanded=True):
        cfg.dea_groupby = st.selectbox("Compare by", ["condition", "sample"], 0)
        cfg.min_cells_per_group = st.number_input("Min cells / group", 5, 1000, 30)
        cfg.dea_logfc_min = st.number_input("Volcano |log2FC| threshold", 0.0, 5.0, 0.25, 0.05)

with pc3:
    with st.expander("Spatial (Squidpy)", expanded=True):
        cfg.spatial_n_neighs = st.slider("Spatial neighbors", 3, 20, 6)
        cfg.cooccurrence_interval = st.slider("Co-occurrence bins", 10, 100, 50)
        cfg.moran_n_perms = st.slider("Moran permutations", 0, 500, 100, 10)
    with st.expander("Report", expanded=False):
        cfg.report_title = st.text_input("Report title", "MACSima Analyzer - Interactive Report")
        cfg.max_cells_in_scatter = st.number_input("Max cells in scatter", 5000, 200000, 60000, 5000)

# --------------------------------------------------------------------------
# Workflow stepper + run
# --------------------------------------------------------------------------
st.subheader("5 | Run workflow")
stepper = st.container()


def draw_stepper(status):
    cols = stepper.columns(len(pipeline.STAGES))
    for col, (sid, label) in zip(cols, pipeline.STAGES):
        s = status.get(sid, {})
        if s.get("done"):
            col.success(f"{label}\n\n{s.get('msg','')}\n\n{s.get('dt',0):.1f}s")
        elif s.get("running"):
            col.info(f"{label}")
        else:
            col.write(f"{label}")


draw_stepper(ss.stage_status)

run = st.button("Run pipeline", type="primary", use_container_width=True,
                disabled=ss.loaded_adata is None)

if run:
    ss.stage_status = {sid: {} for sid, _ in pipeline.STAGES}
    prog = st.progress(0.0, text="Starting...")
    done_count = [0]
    n_stages = len(pipeline.STAGES)

    def cb(stage, msg, dt):
        ss.stage_status[stage] = {"done": True, "msg": msg, "dt": dt}
        done_count[0] += 1
        prog.progress(done_count[0] / n_stages, text=f"{STAGE_LABELS.get(stage, stage)} - {msg}")

    try:
        adata, results = pipeline.run_pipeline(
            None, cfg, condition_map=cond_map, adata=ss.loaded_adata.copy(), progress=cb)
        ss.adata, ss.results = adata, results
        rp = os.path.join(tempfile.gettempdir(), "MACSima_report.html")
        report.build_report(adata, results, cfg, rp)
        ss.report_path = rp
        prog.progress(1.0, text="Done.")
        st.success("Pipeline complete.")
    except Exception as e:
        st.exception(e)
    draw_stepper(ss.stage_status)

# --------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------
if ss.results is not None:
    st.subheader("6 | Results")
    adata, res = ss.adata, ss.results
    r1, r2 = st.columns(2)
    if "cell_type_auto" in adata.obs:
        r1.write("**Auto cell-type counts**")
        r1.dataframe(adata.obs["cell_type_auto"].value_counts(), use_container_width=True)
    dea = res.get("dea")
    if dea and dea.get("single_cell"):
        disease = next((g for g in dea["single_cell"] if "myeloma" in g.lower()),
                       list(dea["single_cell"])[0])
        r2.write(f"**Top up-regulated in {disease} (single-cell)**")
        r2.dataframe(dea["single_cell"][disease].head(10)[
            ["marker", "log2FC", "padj"]], use_container_width=True)

    if ss.report_path and os.path.exists(ss.report_path):
        with open(ss.report_path, "rb") as f:
            st.download_button("Download interactive HTML report", f,
                               file_name="MACSima_report.html", mime="text/html",
                               type="primary", use_container_width=True)
        with st.expander("Preview report inline", expanded=False):
            st.components.v1.html(open(ss.report_path).read(), height=800, scrolling=True)
