#!/usr/bin/env python3
"""
Generate GC-content stratification BED files for a chromosome FASTA.

Usage:
    python generate_gc_strata.py <ref_fasta> <output_dir> [window_size]

Example:
    python evaluation/generate_gc_strata.py \
        data/reference/chr22.fa \
        data/reference/gc_strata \
        1000

Outputs BED files per GC bin and updates stratification_chr22.tsv.
"""

import os
import sys
from pathlib import Path


def read_fasta_simple(fasta_path: str) -> dict:
    """Read single-sequence FASTA without external dependencies."""
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_name is not None:
        sequences[current_name] = "".join(current_seq)

    return sequences


def compute_gc_windows(sequence: str, chrom: str, window_size: int) -> list:
    """Compute GC% per non-overlapping window."""
    windows = []
    seq_len = len(sequence)

    for start in range(0, seq_len, window_size):
        end = min(start + window_size, seq_len)
        subseq = sequence[start:end]

        # Skip windows that are mostly N
        n_count = subseq.count("N")
        valid_bases = len(subseq) - n_count
        if valid_bases < window_size * 0.5:
            continue

        gc_count = subseq.count("G") + subseq.count("C")
        gc_pct = gc_count / valid_bases * 100.0

        windows.append((chrom, start, end, gc_pct))

    return windows


# GC bins matching the reference repo's pattern
GC_BINS = [
    ("GC_0_20",   0,  20),
    ("GC_20_30", 20,  30),
    ("GC_30_40", 30,  40),
    ("GC_40_50", 40,  50),
    ("GC_50_60", 50,  60),
    ("GC_60_80", 60,  80),
    ("GC_80_100", 80, 100),
]


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <ref_fasta> <output_dir> [window_size]")
        sys.exit(1)

    fasta_path = sys.argv[1]
    output_dir = sys.argv[2]
    window_size = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

    os.makedirs(output_dir, exist_ok=True)

    print(f"Reading: {fasta_path}")
    sequences = read_fasta_simple(fasta_path)

    if not sequences:
        print("ERROR: No sequences found in FASTA.")
        sys.exit(1)

    chrom = list(sequences.keys())[0]
    seq = sequences[chrom]
    print(f"  Chromosome: {chrom} ({len(seq):,} bp)")
    print(f"  Window size: {window_size} bp")

    windows = compute_gc_windows(seq, chrom, window_size)
    print(f"  Valid windows: {len(windows)}")

    # Write BED files per bin
    bed_files = {}
    total_written = 0
    for bin_name, gc_min, gc_max in GC_BINS:
        bed_path = os.path.join(output_dir, f"{chrom}_{bin_name}.bed")
        count = 0
        with open(bed_path, "w") as f:
            for c, start, end, gc_pct in windows:
                if gc_min <= gc_pct < gc_max:
                    f.write(f"{c}\t{start}\t{end}\n")
                    count += 1
        bed_files[bin_name] = bed_path
        total_written += count
        print(f"  {bin_name}: {count} windows -> {bed_path}")

    # Write summary
    print(f"\n  Total windows written: {total_written}")

    # Generate/update stratification TSV
    ref_dir = str(Path(fasta_path).parent)
    strat_tsv = os.path.join(ref_dir, "stratification_chr22.tsv")

    # Read existing entries
    existing_entries = []
    if os.path.exists(strat_tsv):
        with open(strat_tsv) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) >= 2 and not parts[0].startswith("GC_"):
                        existing_entries.append(line)

    # Write updated TSV
    with open(strat_tsv, "w") as f:
        for entry in existing_entries:
            f.write(entry + "\n")
        for bin_name, _, _ in GC_BINS:
            bed_rel = f"gc_strata/{chrom}_{bin_name}.bed"
            f.write(f"{bin_name}\t{bed_rel}\n")

    print(f"\n  Updated stratification TSV: {strat_tsv}")
    print("  Done!")


if __name__ == "__main__":
    main()
