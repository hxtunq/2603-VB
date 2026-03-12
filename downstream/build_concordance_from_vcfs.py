#!/usr/bin/env python3
"""Build variant-level concordance files from tp/fp/fn VCFs.

This script discovers caller-specific vcfeval outputs and generates:
- variant_concordance_<cov>x.tsv
- upset_long_callset_<cov>x.csv (FP long format)
- upset_long_baseline_<cov>x.csv (FN long format)

Use this when concordance files are missing but tp/fp/fn.vcf.gz exist.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd


CALLER_ALIAS = {
    "gatk": "HC",
    "hc": "HC",
    "haplotypecaller": "HC",
    "deepvariant": "DV",
    "dv": "DV",
    "strelka2": "ST",
    "strelka": "ST",
    "st": "ST",
    "freebayes": "FB",
    "fb": "FB",
    "dnascope": "DS",
    "ds": "DS",
}

PREFERRED_ALIAS_ORDER = ["HC", "DV", "ST", "FB", "DS"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build concordance tables from tp/fp/fn VCFs")
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument(
        "--search-root",
        type=Path,
        default=Path("results/eval"),
        help="Root to recursively search for tp/fp/fn.vcf.gz",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/eval/vcfeval"),
        help="Directory to write variant_concordance_* and upset_long_*",
    )
    p.add_argument("--coverages", default="10,20,30,50")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def parse_coverages(raw: str) -> List[int]:
    vals = []
    for x in raw.split(","):
        x = x.strip().lower().replace("x", "")
        if x:
            vals.append(int(x))
    return sorted(set(vals))


def find_coverage_in_path(path: Path, allowed: Set[str]) -> str:
    for part in path.parts:
        m = re.fullmatch(r"(\d+)x", part.lower())
        if m and part.lower() in allowed:
            return part.lower()
    return ""


def infer_alias_from_path(path: Path) -> str:
    candidates = [path.parent.name, path.parent.parent.name if path.parent.parent else ""]
    for raw in candidates:
        key = raw.strip().lower()
        if key in CALLER_ALIAS:
            return CALLER_ALIAS[key]
    return ""


def extract_variant_ids(vcf_gz: Path) -> Set[str]:
    out: Set[str] = set()
    with gzip.open(vcf_gz, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            out.add(f"{cols[0]}:{cols[1]}:{cols[3]}:{cols[4]}")
    return out


def discover_metric_files(search_root: Path, coverages: List[int]) -> Dict[Tuple[str, str, str], Path]:
    wanted_cov = {f"{c}x" for c in coverages}
    found: Dict[Tuple[str, str, str], Path] = {}

    for metric in ["tp", "fp", "fn"]:
        for path in search_root.rglob(f"{metric}.vcf.gz"):
            cov = find_coverage_in_path(path, wanted_cov)
            if not cov:
                continue
            alias = infer_alias_from_path(path)
            if not alias:
                continue

            key = (cov, alias, metric)
            if key in found:
                # Prefer paths that include vcfeval token.
                old = str(found[key]).lower()
                new = str(path).lower()
                if "vcfeval" in new and "vcfeval" not in old:
                    found[key] = path
            else:
                found[key] = path

    return found


def build_one_coverage(cov: str, metric_paths: Dict[Tuple[str, str, str], Path], out_dir: Path, overwrite: bool) -> None:
    aliases = sorted({a for (c, a, _m) in metric_paths if c == cov}, key=lambda x: PREFERRED_ALIAS_ORDER.index(x) if x in PREFERRED_ALIAS_ORDER else 999)
    if len(aliases) < 2:
        print(f"[WARN] {cov}: found fewer than 2 callers with tp/fp/fn")
        return

    out_conc = out_dir / f"variant_concordance_{cov}.tsv"
    out_fp = out_dir / f"upset_long_callset_{cov}.csv"
    out_fn = out_dir / f"upset_long_baseline_{cov}.csv"

    if (not overwrite) and out_conc.exists() and out_fp.exists() and out_fn.exists():
        print(f"[SKIP] {cov}: outputs already exist")
        return

    tp: Dict[str, Set[str]] = {}
    fp: Dict[str, Set[str]] = {}
    fn: Dict[str, Set[str]] = {}

    for alias in aliases:
        tp_p = metric_paths.get((cov, alias, "tp"))
        fp_p = metric_paths.get((cov, alias, "fp"))
        fn_p = metric_paths.get((cov, alias, "fn"))
        if not (tp_p and fp_p and fn_p):
            print(f"[WARN] {cov} {alias}: missing one of tp/fp/fn, caller skipped")
            continue
        tp[alias] = extract_variant_ids(tp_p)
        fp[alias] = extract_variant_ids(fp_p)
        fn[alias] = extract_variant_ids(fn_p)

    aliases = [a for a in aliases if a in tp and a in fp and a in fn]
    if len(aliases) < 2:
        print(f"[WARN] {cov}: fewer than 2 complete callers after filtering")
        return

    all_truth = set().union(*[tp[a] | fn[a] for a in aliases]) if aliases else set()
    all_called = set().union(*[tp[a] | fp[a] for a in aliases]) if aliases else set()
    all_variants = sorted(all_truth | all_called)

    rows = []
    for vid in all_variants:
        row = {"variant_id": vid, "truth": int(vid in all_truth)}
        for a in aliases:
            row[a] = int(vid in (tp[a] | fp[a]))
        rows.append(row)

    df = pd.DataFrame(rows)
    cols = ["variant_id", "truth"] + aliases
    df = df[cols]
    df.to_csv(out_conc, sep="\t", index=False)

    with out_fp.open("w", encoding="utf-8") as fh:
        for a in aliases:
            for vid in sorted(fp[a]):
                fh.write(f"{a},{vid}\n")

    with out_fn.open("w", encoding="utf-8") as fh:
        for a in aliases:
            for vid in sorted(fn[a]):
                fh.write(f"{a},{vid}\n")

    print(f"[OK] {out_conc}")
    print(f"[OK] {out_fp}")
    print(f"[OK] {out_fn}")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    search_root = (root / args.search_root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    coverages = parse_coverages(args.coverages)
    metric_paths = discover_metric_files(search_root, coverages)

    if not metric_paths:
        print(f"[WARN] No tp/fp/fn.vcf.gz discovered under {search_root}")
        return

    print(f"[INFO] Discovered {len(metric_paths)} metric files under {search_root}")
    for cov in [f"{c}x" for c in coverages]:
        build_one_coverage(cov, metric_paths, out_dir, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
