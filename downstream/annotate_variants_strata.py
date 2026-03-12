#!/usr/bin/env python3
"""Annotate extracted error-pattern variants with strata information.

Inputs:
- results/analysis/error_patterns/*.tsv
- data/reference/gc_strata/*.bed (GC and CDS_GC bins)
- optional CDS canonical BED
- data/reference/repeats/*.bed (tandem repeat, mappability, segdup)

Outputs:
- results/analysis/stratification/*_with_strata.tsv
- results/analysis/stratification/strata_counts_summary.tsv
"""

from __future__ import annotations

import argparse
import bisect
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


BED_REC = Tuple[int, int, str]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Annotate downstream error patterns with strata")
    p.add_argument("--project-root", type=Path, default=Path("."))
    p.add_argument("--patterns-dir", type=Path, default=Path("results/analysis/error_patterns"))
    p.add_argument("--gc-dir", type=Path, default=Path("data/reference/gc_strata"))
    p.add_argument(
        "--cds-bed",
        type=Path,
        default=Path("data/reference/stratification/CDS-canonical.chr22.bed"),
        help="CDS canonical BED; if missing, in_cds_canonical will be NA",
    )
    p.add_argument("--out-dir", type=Path, default=Path("results/analysis/stratification"))
    p.add_argument(
        "--repeats-dir",
        type=Path,
        default=Path("data/reference/repeats"),
        help="Directory with GIAB repeat/mappability/segdup BEDs",
    )
    return p.parse_args()


def load_bed(path: Path, label: str) -> Dict[str, List[BED_REC]]:
    out: Dict[str, List[BED_REC]] = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            out.setdefault(chrom, []).append((start, end, label))
    for chrom in out:
        out[chrom].sort(key=lambda x: (x[0], x[1]))
    return out


def merge_interval_maps(a: Dict[str, List[BED_REC]], b: Dict[str, List[BED_REC]]) -> Dict[str, List[BED_REC]]:
    out = {k: list(v) for k, v in a.items()}
    for chrom, recs in b.items():
        out.setdefault(chrom, []).extend(recs)
    for chrom in out:
        out[chrom].sort(key=lambda x: (x[0], x[1]))
    return out


def overlap_labels(intervals: List[BED_REC], pos0: int) -> List[str]:
    if not intervals:
        return []
    starts = [x[0] for x in intervals]
    idx = bisect.bisect_right(starts, pos0)
    labels: List[str] = []

    i = idx - 1
    while i >= 0 and intervals[i][0] <= pos0:
        start, end, label = intervals[i]
        if start <= pos0 < end:
            labels.append(label)
        if pos0 - start > 2000000:
            break
        i -= 1

    return labels


def parse_variant_id(v: str) -> Tuple[str, int, str, str]:
    # Expected: CHROM:POS:REF:ALT
    parts = v.split(":")
    if len(parts) != 4:
        raise ValueError(f"Unexpected variant_id format: {v}")
    chrom, pos_s, ref, alt = parts
    return chrom, int(pos_s), ref, alt


def infer_gc_bin(labels: List[str]) -> str:
    gc_only = [x for x in labels if re.search(r"^GC_", x)]
    if not gc_only:
        return ""
    return sorted(gc_only)[0]


def infer_cds_gc_bin(labels: List[str]) -> str:
    cds_gc = [x for x in labels if re.search(r"^CDS_GC_", x)]
    if not cds_gc:
        return ""
    return sorted(cds_gc)[0]


def to_gc_pct_mid(gc_bin: str) -> float:
    # GC_20_30 -> 25.0
    if not gc_bin:
        return float("nan")
    m = re.match(r"^GC_(\d+)_(\d+)$", gc_bin)
    if not m:
        return float("nan")
    lo = float(m.group(1))
    hi = float(m.group(2))
    return (lo + hi) / 2.0


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    patterns_dir = (root / args.patterns_dir).resolve()
    gc_dir = (root / args.gc_dir).resolve()
    cds_bed = (root / args.cds_bed).resolve()
    repeats_dir = (root / args.repeats_dir).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    gc_map: Dict[str, List[BED_REC]] = {}
    for bed in sorted(gc_dir.glob("*.bed")):
        label = bed.stem.split("_", 1)[1] if "_" in bed.stem else bed.stem
        gc_map = merge_interval_maps(gc_map, load_bed(bed, label))

    cds_map = load_bed(cds_bed, "CDS_CANONICAL")

    # GIAB repeat / mappability / segdup BEDs
    tandem_map = load_bed(repeats_dir / "tandem_repeat_homopolymer.bed", "tandem_repeat")
    lowmapq_map = load_bed(repeats_dir / "low_mappability.bed", "low_mapq")
    segdup_map = load_bed(repeats_dir / "segdup.bed", "segdup")

    files = [
        p for p in sorted(patterns_dir.glob("*.tsv"))
        if re.search(r"_\d+x\.tsv$", p.name)
    ]
    if not files:
        print(f"[WARN] No pattern files in {patterns_dir}")
        return

    summary_rows = []

    for path in files:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            out_file = out_dir / path.name.replace(".tsv", "_with_strata.tsv")
            df.to_csv(out_file, sep="\t", index=False)
            continue

        chroms = []
        pos = []
        refs = []
        alts = []
        gc_bins = []
        cds_gc_bins = []
        in_cds = []
        in_tandem = []
        in_lowmapq = []
        in_segdup = []

        for vid in df["variant_id"].astype(str):
            c, p, r, a = parse_variant_id(vid)
            chroms.append(c)
            pos.append(p)
            refs.append(r)
            alts.append(a)
            pos0 = p - 1

            labels = overlap_labels(gc_map.get(c, []), pos0)
            gc_bin = infer_gc_bin(labels)
            cds_gc_bin = infer_cds_gc_bin(labels)
            gc_bins.append(gc_bin)
            cds_gc_bins.append(cds_gc_bin)

            cds_labels = overlap_labels(cds_map.get(c, []), pos0)
            in_cds.append(1 if cds_labels else 0)

            in_tandem.append(1 if overlap_labels(tandem_map.get(c, []), pos0) else 0)
            in_lowmapq.append(1 if overlap_labels(lowmapq_map.get(c, []), pos0) else 0)
            in_segdup.append(1 if overlap_labels(segdup_map.get(c, []), pos0) else 0)

        out = df.copy()
        out["chrom"] = chroms
        out["pos"] = pos
        out["ref"] = refs
        out["alt"] = alts
        out["gc_bin"] = gc_bins
        out["gc_pct_mid"] = [to_gc_pct_mid(x) for x in gc_bins]
        out["cds_gc_bin"] = cds_gc_bins
        out["in_cds_canonical"] = in_cds
        # Protein coding fallback: use CDS canonical when no separate BED is provided.
        out["in_protein_coding"] = in_cds
        out["in_tandem_repeat"] = in_tandem
        out["in_low_mapq"] = in_lowmapq
        out["in_segdup"] = in_segdup

        out_file = out_dir / path.name.replace(".tsv", "_with_strata.tsv")
        out.to_csv(out_file, sep="\t", index=False)
        print(f"[OK] {out_file}")

        patt = out["pattern"].iloc[0] if "pattern" in out.columns and len(out) else ""
        cov_match = re.search(r"_(\d+x)\.tsv$", path.name)
        coverage = cov_match.group(1) if cov_match else ""
        summary_rows.append(
            {
                "pattern": patt,
                "coverage": coverage,
                "n_variants": int(len(out)),
                "n_in_cds": int(out["in_cds_canonical"].sum()),
                "n_in_tandem_repeat": int(out["in_tandem_repeat"].sum()),
                "n_in_low_mapq": int(out["in_low_mapq"].sum()),
                "n_in_segdup": int(out["in_segdup"].sum()),
                "median_gc_pct_mid": float(out["gc_pct_mid"].median(skipna=True)),
            }
        )

    if summary_rows:
        summary = pd.DataFrame(summary_rows)
        summary_file = out_dir / "strata_counts_summary.tsv"
        summary.to_csv(summary_file, sep="\t", index=False)
        print(f"[OK] Summary -> {summary_file}")


if __name__ == "__main__":
    main()
