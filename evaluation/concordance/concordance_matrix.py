#!/usr/bin/env python3
"""
Build a concordance matrix from pairwise RTG vcfeval results.

Usage:
    python concordance_matrix.py <concordance_root> <output_tsv>

Example:
    python evaluation/concordance/concordance_matrix.py \
        results/eval/concordance \
        results/eval/concordance/concordance_matrix.tsv

The script reads summary.txt from each pair directory and computes:
    concordance_rate = TP / (TP + FP + FN)
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path


CALLER_ORDER = ["HC", "DV", "ST", "FB", "DS"]


def parse_summary(summary_path: str) -> dict:
    """Parse RTG vcfeval summary.txt and return TP, FP, FN counts."""
    results = {}
    with open(summary_path) as f:
        lines = f.readlines()

    # Find header line
    header_idx = None
    for i, line in enumerate(lines):
        if "Threshold" in line:
            header_idx = i
            break

    if header_idx is None:
        return results

    # Parse the "None" threshold line (unfiltered) — usually the last data line
    for line in lines[header_idx + 1:]:
        line = line.strip()
        if not line or line.startswith("-"):
            continue
        parts = line.split()
        if len(parts) >= 7:
            results = {
                "threshold": parts[0],
                "tp_baseline": float(parts[1]),
                "tp_call": float(parts[2]),
                "fp": float(parts[3]),
                "fn": float(parts[4]),
                "precision": float(parts[5]) if parts[5] != "N/A" else 0.0,
                "recall": float(parts[6]) if parts[6] != "N/A" else 0.0,
            }

    return results


def compute_concordance(tp: float, fp: float, fn: float) -> float:
    """Concordance rate = TP / (TP + FP + FN)."""
    denom = tp + fp + fn
    if denom == 0:
        return 0.0
    return tp / denom


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <concordance_root> <output_tsv>")
        sys.exit(1)

    concord_root = Path(sys.argv[1])
    output_tsv = sys.argv[2]

    # Collect results:  coverage -> (callerA, callerB) -> concordance
    all_results = defaultdict(dict)

    for cov_dir in sorted(concord_root.iterdir()):
        if not cov_dir.is_dir() or not re.match(r"\d+x$", cov_dir.name):
            continue
        cov = cov_dir.name

        for pair_dir in sorted(cov_dir.iterdir()):
            if not pair_dir.is_dir():
                continue
            match = re.match(r"^([A-Z]+)_vs_([A-Z]+)$", pair_dir.name)
            if not match:
                continue
            caller_a, caller_b = match.group(1), match.group(2)

            summary = pair_dir / "summary.txt"
            if not summary.exists():
                print(f"  WARN: No summary.txt in {pair_dir}")
                continue

            stats = parse_summary(str(summary))
            if not stats:
                print(f"  WARN: Could not parse {summary}")
                continue

            tp = stats["tp_baseline"]
            fp = stats["fp"]
            fn = stats["fn"]
            rate = compute_concordance(tp, fp, fn)

            all_results[cov][(caller_a, caller_b)] = rate
            all_results[cov][(caller_b, caller_a)] = rate  # symmetric

    # Write long-form TSV
    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    with open(output_tsv, "w") as f:
        f.write("Coverage\tCallerA\tCallerB\tConcordance\n")
        for cov in sorted(all_results.keys(), key=lambda x: int(x.replace("x", ""))):
            for (a, b), rate in sorted(all_results[cov].items()):
                f.write(f"{cov}\t{a}\t{b}\t{rate:.6f}\n")
            # Add self-concordance = 1.0
            for caller in CALLER_ORDER:
                f.write(f"{cov}\t{caller}\t{caller}\t1.000000\n")

    print(f"Concordance matrix written to: {output_tsv}")
    print(f"Coverages: {sorted(all_results.keys())}")

    # Print a summary per coverage
    for cov in sorted(all_results.keys(), key=lambda x: int(x.replace("x", ""))):
        print(f"\n--- {cov} ---")
        pairs = all_results[cov]
        seen = set()
        for (a, b), rate in sorted(pairs.items()):
            key = tuple(sorted([a, b]))
            if key not in seen:
                seen.add(key)
                print(f"  {a} vs {b}: {rate:.4f}")


if __name__ == "__main__":
    main()
