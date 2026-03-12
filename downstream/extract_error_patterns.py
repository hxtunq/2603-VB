#!/usr/bin/env python3
"""Extract FN/FP intersection patterns from variant concordance tables.

This script does NOT modify existing pipeline scripts. It reads:
  results/eval/vcfeval/variant_concordance_{cov}x.tsv

And writes downstream pattern tables:
  results/analysis/error_patterns/<pattern>_{cov}x.tsv

Patterns:
- all_5_fn: truth=1 and all 5 callers missed
- all_5_fp: truth=0 and all 5 callers called FP
- single_caller_fn: truth=1 and exactly 1 caller missed
- single_caller_fp: truth=0 and exactly 1 caller called FP
- dv_ds_fn: truth=1 and both DV and DS missed (may overlap all_5_fn)
- dv_ds_fp: truth=0 and both DV and DS called FP (may overlap all_5_fp)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


DEFAULT_ALIAS = ["HC", "DV", "ST", "FB", "DS"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract FN/FP error patterns from concordance tables")
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument(
        "--vcfeval-dir",
        type=Path,
        default=Path("results/eval/vcfeval"),
        help="Directory containing variant_concordance_<cov>x.tsv files",
    )
    p.add_argument("--out-dir", type=Path, default=Path("results/analysis/error_patterns"))
    p.add_argument("--coverages", default="10,20,30,50", help="Comma-separated coverages, e.g. 10,20,30,50")
    p.add_argument("--callers", default=",".join(DEFAULT_ALIAS), help="Caller alias columns in concordance table")
    return p.parse_args()


def parse_coverages(raw: str) -> List[int]:
    vals = []
    for x in raw.split(","):
        x = x.strip().lower().replace("x", "")
        if not x:
            continue
        vals.append(int(x))
    return sorted(set(vals))


def safe_focus_caller(row: pd.Series, caller_cols: List[str], mode: str) -> str:
    if mode == "fn":
        missed = [c for c in caller_cols if int(row[c]) == 0]
        return missed[0] if len(missed) == 1 else ""
    called = [c for c in caller_cols if int(row[c]) == 1]
    return called[0] if len(called) == 1 else ""


def extract_for_coverage(df: pd.DataFrame, caller_cols: List[str]) -> Dict[str, pd.DataFrame]:
    work = df.copy()
    work["n_called"] = work[caller_cols].sum(axis=1)
    n_callers = len(caller_cols)

    all_5_fn = work[(work["truth"] == 1) & (work["n_called"] == 0)].copy()
    all_5_fn["pattern"] = "all_5_fn"
    all_5_fn["focus_caller"] = ""

    all_5_fp = work[(work["truth"] == 0) & (work["n_called"] == n_callers)].copy()
    all_5_fp["pattern"] = "all_5_fp"
    all_5_fp["focus_caller"] = ""

    single_fn = work[(work["truth"] == 1) & (work["n_called"] == n_callers - 1)].copy()
    single_fn["pattern"] = "single_caller_fn"
    single_fn["focus_caller"] = single_fn.apply(
        lambda r: safe_focus_caller(r, caller_cols, "fn"), axis=1
    )

    single_fp = work[(work["truth"] == 0) & (work["n_called"] == 1)].copy()
    single_fp["pattern"] = "single_caller_fp"
    single_fp["focus_caller"] = single_fp.apply(
        lambda r: safe_focus_caller(r, caller_cols, "fp"), axis=1
    )

    if "DV" in caller_cols and "DS" in caller_cols:
        dv_ds_fn = work[(work["truth"] == 1) & (work["DV"] == 0) & (work["DS"] == 0)].copy()
        dv_ds_fn["pattern"] = "dv_ds_fn"
        dv_ds_fn["focus_caller"] = "DV,DS"
        dv_ds_fn["exclusive_to_dv_ds"] = (dv_ds_fn["n_called"] == n_callers - 2).astype(int)

        dv_ds_fp = work[(work["truth"] == 0) & (work["DV"] == 1) & (work["DS"] == 1)].copy()
        dv_ds_fp["pattern"] = "dv_ds_fp"
        dv_ds_fp["focus_caller"] = "DV,DS"
        dv_ds_fp["exclusive_to_dv_ds"] = (dv_ds_fp["n_called"] == 2).astype(int)
    else:
        dv_ds_fn = work.iloc[0:0].copy()
        dv_ds_fp = work.iloc[0:0].copy()

    for frame in [all_5_fn, all_5_fp, single_fn, single_fp, dv_ds_fn, dv_ds_fp]:
        if "exclusive_to_dv_ds" not in frame.columns:
            frame["exclusive_to_dv_ds"] = 0

    return {
        "all_5_fn": all_5_fn,
        "all_5_fp": all_5_fp,
        "single_caller_fn": single_fn,
        "single_caller_fp": single_fp,
        "dv_ds_fn": dv_ds_fn,
        "dv_ds_fp": dv_ds_fp,
    }


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    vcfeval_dir = (project_root / args.vcfeval_dir).resolve()
    out_dir = (project_root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    coverages = parse_coverages(args.coverages)
    preferred_callers = [x.strip() for x in args.callers.split(",") if x.strip()]

    summary_rows = []

    for cov in coverages:
        concordance = vcfeval_dir / f"variant_concordance_{cov}x.tsv"
        if not concordance.exists():
            print(f"[WARN] Missing concordance file: {concordance}")
            continue

        df = pd.read_csv(concordance, sep="\t")
        missing = {"variant_id", "truth"} - set(df.columns)
        if missing:
            if {"Coverage", "CallerA", "CallerB", "Concordance"}.issubset(df.columns):
                raise ValueError(
                    "Detected pairwise concordance summary table, not variant-level concordance. "
                    "This pipeline requires files like variant_concordance_10x.tsv with columns: "
                    "variant_id, truth, HC, DV, ST, FB, DS."
                )
            raise ValueError(f"Missing required columns {missing} in {concordance}")

        caller_cols = [c for c in preferred_callers if c in df.columns]
        if len(caller_cols) < 2:
            raise ValueError(
                f"Not enough caller columns in {concordance}. Found {caller_cols}"
            )

        patterns = extract_for_coverage(df, caller_cols)

        base_cols = ["variant_id", "truth"] + caller_cols + [
            "n_called", "pattern", "focus_caller", "exclusive_to_dv_ds"
        ]

        for pattern_name, frame in patterns.items():
            out_file = out_dir / f"{pattern_name}_{cov}x.tsv"
            frame = frame[base_cols].sort_values(["truth", "n_called", "variant_id"])
            frame.to_csv(out_file, sep="\t", index=False)
            summary_rows.append(
                {
                    "coverage": f"{cov}x",
                    "pattern": pattern_name,
                    "n_variants": int(len(frame)),
                }
            )
            print(f"[OK] {out_file} ({len(frame)} variants)")

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows).sort_values(["coverage", "pattern"])
        summary_file = out_dir / "pattern_counts_summary.tsv"
        summary_df.to_csv(summary_file, sep="\t", index=False)
        print(f"[OK] Summary -> {summary_file}")
    else:
        print("[WARN] No outputs generated. Check input concordance files.")


if __name__ == "__main__":
    main()
