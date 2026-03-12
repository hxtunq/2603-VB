#!/usr/bin/env python3
"""
Generate GC-content stratification BED files specifically for CDS regions.

Usage:
    python generate_cds_gc_strata.py <ref_fasta> <cds_bed> <output_dir>

Outputs BED files per GC bin for the CDS regions and updates stratification_chr22.tsv.
"""

import os
import sys
from pathlib import Path

def read_fasta_simple(fasta_path: str) -> dict:
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

GC_BINS = [
    ("CDS_GC_0_20",   0,  20),
    ("CDS_GC_20_30", 20,  30),
    ("CDS_GC_30_40", 30,  40),
    ("CDS_GC_40_50", 40,  50),
    ("CDS_GC_50_60", 50,  60),
    ("CDS_GC_60_80", 60,  80),
    ("CDS_GC_80_100", 80, 100),
]

def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <ref_fasta> <cds_bed> <output_dir>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    cds_bed = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    print(f"Reading FASTA: {fasta_path}")
    sequences = read_fasta_simple(fasta_path)

    print(f"Reading BED: {cds_bed}")
    intervals = []
    with open(cds_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            intervals.append((chrom, start, end))
            
    chrom_out = intervals[0][0] if intervals else list(sequences.keys())[0]

    # Categorize
    binned_intervals = {name: [] for name, _, _ in GC_BINS}
    
    for chrom, start, end in intervals:
        if chrom not in sequences:
            continue
        seq = sequences[chrom][start:end]
        n_count = seq.count("N")
        valid_bases = len(seq) - n_count
        if valid_bases == 0:
            continue
            
        gc_count = seq.count("G") + seq.count("C")
        gc_pct = gc_count / valid_bases * 100.0
        
        for name, gc_min, gc_max in GC_BINS:
            if gc_min <= gc_pct < gc_max:
                binned_intervals[name].append((chrom, start, end))
                break

    total_written = 0
    for name in binned_intervals:
        bed_path = os.path.join(output_dir, f"{chrom_out}_{name}.bed")
        count = len(binned_intervals[name])
        with open(bed_path, "w") as f:
            for chrom, start, end in binned_intervals[name]:
                f.write(f"{chrom}\t{start}\t{end}\n")
        total_written += count
        print(f"  {name}: {count} intervals -> {bed_path}")

    # Update TSV
    ref_dir = str(Path(fasta_path).parent)
    strat_tsv = os.path.join(ref_dir, "stratification_chr22.tsv")

    existing_entries = []
    if os.path.exists(strat_tsv):
        with open(strat_tsv) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) >= 2 and not parts[0].startswith("CDS_GC_"):
                        existing_entries.append(line)

    with open(strat_tsv, "w") as f:
        for entry in existing_entries:
            f.write(entry + "\n")
        for name, _, _ in GC_BINS:
            bed_rel = f"gc_strata/{chrom_out}_{name}.bed"
            f.write(f"{name}\t{bed_rel}\n")

    print(f"\n  Updated stratification TSV: {strat_tsv}")
    print("  Done!")

if __name__ == "__main__":
    main()
