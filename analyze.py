"""Headless runner - analyse a folder of MACSima CSV exports and build the report.

This is the software's command-line entry point; it operates on real data only.
(For the synthetic benchmark, see benchmark/run_benchmark.py.)

Examples
--------
python analyze.py --data /path/to/csv_folder --out report.html
python analyze.py --data ./my_csvs --out report.html \
    --cofactor 5 --resolution 1.0 --groupby condition
"""
from __future__ import annotations
import argparse
import os
import time

from macsima.config import PipelineConfig
from macsima import pipeline, report


def main():
    ap = argparse.ArgumentParser(description="MACSima Analyzer - headless analysis")
    ap.add_argument("--data", required=True,
                    help="folder containing MACSima CSV files (one per sample)")
    ap.add_argument("--out", default="report.html", help="output HTML report path")
    ap.add_argument("--cofactor", type=float, default=5.0)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--groupby", default="condition", choices=["condition", "sample"])
    ap.add_argument("--no-annotate", action="store_true")
    args = ap.parse_args()

    if not os.path.isdir(args.data):
        ap.error(f"--data must be a folder of CSVs (got: {args.data})")

    cfg = PipelineConfig()
    cfg.cofactor = args.cofactor
    cfg.leiden_resolution = args.resolution
    cfg.n_neighbors = args.n_neighbors
    cfg.dea_groupby = args.groupby
    cfg.auto_annotate = not args.no_annotate

    def prog(stage, msg, dt):
        print(f"  [{stage:9s}] {msg}  ({dt:.1f}s)")

    print(f"Analysing: {args.data}")
    t0 = time.time()
    adata, results = pipeline.run_pipeline(args.data, cfg, progress=prog)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    report.build_report(adata, results, cfg, args.out)
    print(f"\nDone in {time.time()-t0:.1f}s. Report: {args.out}  "
          f"({os.path.getsize(args.out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
