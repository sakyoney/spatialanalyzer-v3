# MACSima Analyzer

Spatial-proteomics analysis for Miltenyi MACSima / MACS iQ View output. It takes
the per-cell object table exported after segmentation, runs a single-cell and
spatial workflow with Scanpy and Squidpy, and produces a self-contained
interactive HTML report. A Streamlit GUI exposes the parameters; the same
pipeline runs from the command line.

The tool works on segmented cell tables (one row per cell, marker mean
intensities plus centroid coordinates), so it sits downstream of image
processing.

## Features

- CSV import with automatic detection of coordinates, area, sample, gate and
  marker columns.
- Channel classification that separates biological markers from DAPI cycle
  channels, bare-fluorophore channels and buffer-only controls.
- arcsinh (or log1p) transform and z-scoring across all samples.
- Per-cell QC and per-marker signal-to-noise against control channels.
- Leiden clustering with UMAP.
- Two cell-type labellings shown together: the gated QiGate labels from the
  export, and unsupervised clusters labelled from a marker-signature dictionary.
- Differential expression by condition, with both a single-cell Wilcoxon test
  and a pseudobulk (per-sample) test, plus differential abundance.
- Squidpy spatial analysis: neighbourhood enrichment, co-occurrence and Moran's
  I, computed on per-sample graphs.
- One interactive HTML report containing all of the above.

## Installation

Python 3.10-3.12.

```bash
pip install -r requirements.txt
```

The main dependencies are scanpy, squidpy, leidenalg/igraph and streamlit.
Static PNG export of figures (optional) needs `kaleido` and Google Chrome; the
HTML report does not.

## Usage

GUI:

```bash
streamlit run app.py
```

Upload your CSVs or point to a folder, review the channel classification, map
samples to conditions, adjust parameters if needed, and run. The report can be
downloaded or previewed inline.

Command line:

```bash
python analyze.py --data /path/to/csv_folder --out report.html
python analyze.py --data ./my_csvs --out report.html --cofactor 5 --resolution 1.0
```

Python:

```python
from macsima.config import PipelineConfig
from macsima import pipeline, report

cfg = PipelineConfig()
adata, results = pipeline.run_pipeline("/path/to/csv_folder", cfg)
report.build_report(adata, results, cfg, "report.html")
```

## Input format

One CSV per sample (or a single CSV with a sample column), one row per cell.
Column names are matched case-insensitively.

| Field | Accepted names | Required |
|---|---|---|
| X centroid | `Centroid X`, `X`, `Position X`, `x_centroid` | yes |
| Y centroid | `Centroid Y`, `Y`, `Position Y`, `y_centroid` | yes |
| Area | any column containing `area` | no |
| Cell id | `Object ID`, `Cell ID`, `ID` | no |
| Sample | `Sample`, `ROI`, `Region`, `Slide` | no (falls back to file name) |
| Condition | `Condition`, `Group`, `Diagnosis`, `Status` | no (set in GUI) |
| Gate | `QiGate`, `Gate`, `Cell type`, `Population` | no |
| Markers | every remaining numeric column | - |

Marker columns may carry a statistic suffix (`CD138 REA929 Mean`), which is
stripped to `CD138 REA929`. Columns whose name starts with an underscore are
carried into `adata.obs` as metadata but not used in analysis.

Channel classes are configurable in `macsima/config.py`, and can be overridden
per run in the GUI.

## Notes on the analysis

- z-scoring is done jointly across all samples after concatenation. Per-sample
  scaling removes between-condition differences and weakens differential
  expression.
- Fold-changes are computed from linear raw intensities so their sign and size
  stay interpretable.
- The single-cell Wilcoxon test pools all cells and has high power but treats
  cells as independent. The pseudobulk test averages each marker per sample and
  compares across samples, which is the appropriate unit for a between-sample
  comparison but is limited by sample count (with 3 vs 3 the smallest two-sided
  Mann-Whitney p is 0.1).
- DAPI and control channels are excluded from clustering and DEA but kept for
  the QC signal-to-noise baseline.
- Spatial graphs use `library_key="sample"` so edges stay within a section.

## Cell-type signatures

Auto-annotation uses the `SIGNATURES` dictionary in `macsima/annotation.py`:

```python
SIGNATURES = {
    "Plasma cell / Myeloma": {"CD138": 1.0, "BCMA": 1.0, "CD38": 0.8,
                              "GPRC5D": 0.9, "CD319": 0.7, "CD45": -0.5},
    "T cell":                {"CD3": 1.0, "CD45": 0.6, "CD20": -0.5},
    # ...
}
```

Keys are substrings matched against channel names (`CD138` matches
`CD138 REA929`); weights may be negative. Each cluster is scored as the weighted
mean z-score of its signature markers and gets the best-scoring label. Edit this
dictionary for a different panel.

## Benchmarking

Synthetic-data generation and validation live under `benchmark/` and are
separate from the software (the software imports nothing from `benchmark/`).

```bash
python benchmark/run_benchmark.py
```

This generates a synthetic Multiple Myeloma dataset with known ground truth,
runs the pipeline on it, scores how well the seeded biology is recovered
(clustering agreement, plasma-cell detection, DEA and spatial marker ranking),
and writes an example report to `benchmark/example_report.html`.

To only generate data:

```bash
python benchmark/generate_mm.py            # -> benchmark/data/synthetic_mm/
```

On the default 3-vs-3 synthetic dataset the pipeline recovers the seeded cell
types and ranks BCMA, GPRC5D, CD138 and CD38 at the top of both the
differential-expression and Moran's I results. Because the data is synthetic
these numbers reflect a controlled sanity check, not performance on real tissue.

## Layout

```
macsima_analyzer/
├── app.py                 # Streamlit GUI
├── analyze.py             # command-line runner
├── requirements.txt
├── macsima/               # analysis package
│   ├── config.py
│   ├── io.py
│   ├── preprocessing.py
│   ├── qc.py
│   ├── annotation.py
│   ├── dea.py
│   ├── spatial.py
│   ├── report.py
│   └── pipeline.py
└── benchmark/             # synthetic data + validation (separate from software)
    ├── generate_mm.py
    └── run_benchmark.py
```

## Disclaimer

This is a research tool. It can highlight patterns for review but does not make
diagnoses and is not a medical device. Clinical interpretation is the
responsibility of a qualified clinician.
