#!/usr/bin/env python3
"""Aggregate AlphaGenome scores onto downstream pattern tables.

Inputs:
- results/analysis/stratification/*_with_strata.tsv
- results/alphagenome/*_variant_scores_tidy.csv

Outputs:
- results/analysis/alphagenome/*_scored.tsv
- results/analysis/alphagenome/scored_summary.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


CALLER_ALIAS_TO_CANONICAL = {
    "HC": "gatk",
    "DV": "deepvariant",
    "ST": "strelka2",
    "FB": "freebayes",
    "DS": "dnascope",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Merge AlphaGenome scores with downstream pattern tables")
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument("--strata-dir", type=Path, default=Path("results/analysis/stratification"))
    p.add_argument("--alphagenome-dir", type=Path, default=Path("results/alphagenome"))
    p.add_argument("--out-dir", type=Path, default=Path("results/analysis/alphagenome"))
    return p.parse_args()


def normalize_variant_id(v: str) -> str:
    # Convert chr22:123:A>T -> chr22:123:A:T when needed.
    if ">" in v and v.count(":") == 2:
        chrom, pos, ra = v.split(":")
        ref, alt = ra.split(">")
        return f"{chrom}:{pos}:{ref}:{alt}"
    return v


def find_numeric_score_cols(df: pd.DataFrame) -> List[str]:
    numeric_cols = []
    blacklist = {"coverage"}
    for c in df.columns:
        if c in blacklist:
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            numeric_cols.append(c)
    return numeric_cols


def build_score_table(alphagenome_dir: Path) -> pd.DataFrame:
    rows = []
    files = sorted(alphagenome_dir.glob("*_variant_scores_tidy.csv"))
    if not files:
        return pd.DataFrame(columns=["variant_id", "caller", "coverage", "error_type", "ag_score"])

    for f in files:
        # filename format: <caller>_<cov>x_<FN|FP>_variant_scores_tidy.csv
        m = re.match(r"^(.+?)_(\d+x)_(FN|FP)_variant_scores_tidy\.csv$", f.name)
        if not m:
            continue
        caller = m.group(1)
        coverage = m.group(2)
        error_type = m.group(3)

        df = pd.read_csv(f)
        if df.empty:
            continue

        if "variant_id" not in df.columns:
            continue

        df["variant_id"] = df["variant_id"].astype(str).map(normalize_variant_id)
        score_cols = find_numeric_score_cols(df)
        if not score_cols:
            continue

        # Robust single score: max absolute numeric signal across score columns.
        arr = df[score_cols].to_numpy(dtype=float)
        score = np.nanmax(np.abs(arr), axis=1)
        tmp = pd.DataFrame(
            {
                "variant_id": df["variant_id"].astype(str),
                "caller": caller,
                "coverage": coverage,
                "error_type": error_type,
                "ag_score": score,
            }
        )

        # Collapse duplicated rows from tidy output by median.
        tmp = (
            tmp.groupby(["variant_id", "caller", "coverage", "error_type"], as_index=False)["ag_score"]
            .median()
        )
        rows.append(tmp)

    if not rows:
        return pd.DataFrame(columns=["variant_id", "caller", "coverage", "error_type", "ag_score"])

    return pd.concat(rows, ignore_index=True)


def expected_callers_for_row(row: pd.Series) -> List[str]:
    patt = str(row.get("pattern", ""))
    focus = str(row.get("focus_caller", ""))

    if patt.startswith("single_caller_") and focus:
        return [CALLER_ALIAS_TO_CANONICAL.get(focus, focus)]
    if patt.startswith("dv_ds_"):
        return ["deepvariant", "dnascope"]

    alias_cols = [a for a in ["HC", "DV", "ST", "FB", "DS"] if a in row.index]
    return [CALLER_ALIAS_TO_CANONICAL[a] for a in alias_cols]


def add_percentile(scores: pd.Series) -> pd.Series:
    if scores.empty:
        return scores
    return scores.rank(pct=True, method="average")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    strata_dir = (root / args.strata_dir).resolve()
    ag_dir = (root / args.alphagenome_dir).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    score_table = build_score_table(ag_dir)
    if score_table.empty:
        print(f"[WARN] No AlphaGenome score files found in {ag_dir}")

    strata_files = sorted(strata_dir.glob("*_with_strata.tsv"))
    if not strata_files:
        print(f"[WARN] No strata files found in {strata_dir}")
        return

    summary_rows: List[Dict[str, object]] = []

    for sf in strata_files:
        m = re.search(r"_(\d+x)_with_strata\.tsv$", sf.name)
        if not m:
            continue
        coverage = m.group(1)

        df = pd.read_csv(sf, sep="\t")
        if df.empty:
            out_file = out_dir / sf.name.replace("_with_strata.tsv", "_scored.tsv")
            df.to_csv(out_file, sep="\t", index=False)
            continue

        if "pattern" not in df.columns:
            raise ValueError(f"Missing pattern column in {sf}")

        df["coverage"] = coverage
        df["error_type"] = np.where(df["pattern"].str.endswith("_fn"), "FN", "FP")

        means = []
        maxs = []
        n_scored = []

        for _, row in df.iterrows():
            callers = expected_callers_for_row(row)
            subset = score_table[
                (score_table["variant_id"] == str(row["variant_id"]))
                & (score_table["coverage"] == coverage)
                & (score_table["error_type"] == str(row["error_type"]))
                & (score_table["caller"].isin(callers))
            ]
            vals = subset["ag_score"].to_numpy(dtype=float)
            if len(vals) == 0:
                means.append(np.nan)
                maxs.append(np.nan)
                n_scored.append(0)
            else:
                means.append(float(np.nanmean(vals)))
                maxs.append(float(np.nanmax(vals)))
                n_scored.append(int(np.sum(~np.isnan(vals))))

        out = df.copy()
        out["alphagenome_score_mean"] = means
        out["alphagenome_score_max"] = maxs
        out["alphagenome_n_callers_scored"] = n_scored
        out["alphagenome_score_percentile"] = add_percentile(out["alphagenome_score_mean"])

        out_file = out_dir / sf.name.replace("_with_strata.tsv", "_scored.tsv")
        out.to_csv(out_file, sep="\t", index=False)
        print(f"[OK] {out_file}")

        summary_rows.append(
            {
                "file": sf.name,
                "n_variants": int(len(out)),
                "n_scored": int((out["alphagenome_n_callers_scored"] > 0).sum()),
                "median_ag_score_mean": float(out["alphagenome_score_mean"].median(skipna=True)),
                "median_ag_score_max": float(out["alphagenome_score_max"].median(skipna=True)),
            }
        )

    if summary_rows:
        summary = pd.DataFrame(summary_rows)
        summary_file = out_dir / "scored_summary.tsv"
        summary.to_csv(summary_file, sep="\t", index=False)
        print(f"[OK] Summary -> {summary_file}")


if __name__ == "__main__":
    main()
