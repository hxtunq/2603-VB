#!/usr/bin/env python3
"""Score FN or FP variants of a SINGLE caller with AlphaGenome.

Each invocation handles exactly one caller × one coverage × one error-type.
Results are written incrementally in batches to avoid OOM.

Supported callers: gatk, deepvariant, strelka2, freebayes, dnascope
Supported coverages: 10, 20, 30, 50 (or 'all' to loop over all)

Inputs:
  results/eval/vcfeval/<COV>x/<caller>/<fn|fp>.vcf.gz

Outputs:
  results/alphagenome/<caller>_<COV>x_<FN|FP>_variant_scores_tidy.csv

Usage:
  export ALPHAGENOME_API_KEY="your-key"

  # Score deepvariant FP at 30x
  python batch_alphagenome_fn_fp.py --caller deepvariant --error-type FP --coverage 30

  # Score gatk FN for all coverages
  python batch_alphagenome_fn_fp.py --caller gatk --error-type FN --coverage all

  # Score dnascope FN at 50x with smaller batches
  python batch_alphagenome_fn_fp.py --caller dnascope --error-type FN --coverage 50 --batch-size 20

  # Custom sequence length / scorers
  python batch_alphagenome_fn_fp.py --caller strelka2 --error-type FP --coverage 30 \
      --sequence-length 500KB --scorers rna_seq,atac,cage
"""

from __future__ import annotations

import argparse
import gc
import gzip
import os
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# ── CLI ──────────────────────────────────────────────────────────────
ALL_COVERAGES = [10, 20, 30, 50]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Score FN/FP variants of a single caller with AlphaGenome"
    )
    p.add_argument("--caller", required=True,
                   help="Caller name: gatk, deepvariant, strelka2, freebayes, dnascope")
    p.add_argument("--error-type", required=True, choices=["FN", "FP"],
                   help="Error type: FN or FP")
    p.add_argument("--coverage", required=True,
                   help="Coverage level (10, 20, 30, 50) or 'all' to loop over all")
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--rtg-dir", type=Path,
                   default=Path("results/eval/vcfeval"))
    p.add_argument("--out-dir", type=Path,
                   default=Path("results/alphagenome"))
    p.add_argument("--organism", default="human",
                   choices=("human", "mouse"))
    p.add_argument("--sequence-length", default="1MB",
                   choices=("16KB", "100KB", "500KB", "1MB"))
    p.add_argument("--scorers", default="all",
                   help="Comma-separated scorer keys, or 'all'")
    p.add_argument("--api-key", default="")
    p.add_argument("--api-key-env", default="ALPHAGENOME_API_KEY")
    p.add_argument("--batch-size", type=int, default=50,
                   help="Number of variants to score per batch (lower = less RAM)")
    p.add_argument("--no-parquet", action="store_true")
    return p.parse_args()


# ── helpers ──────────────────────────────────────────────────────────
def normalize_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def load_vcf_rows(path: Path, caller: str, error_type: str) -> pd.DataFrame:
    """Parse a gzipped VCF into a flat DataFrame of variant metadata."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5 or "," in cols[4]:
                continue
            chrom = normalize_chrom(cols[0])
            pos = int(cols[1])
            ref = cols[3].upper()
            alt = cols[4].upper()
            rows.append({
                "variant_id": f"{chrom}:{pos}:{ref}>{alt}",
                "caller": caller,
                "error_type": error_type,
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
            })
    return pd.DataFrame(rows)


# ── score one caller × one coverage × one error-type ─────────────────
def score_one(args: argparse.Namespace, caller: str, error_type: str,
              cov: int, dna_model, organism, sequence_length,
              selected_scorers, out_dir: Path) -> None:
    """Score FN or FP for a single caller at a specific coverage."""

    root = args.root.resolve()
    rtg_dir = (root / args.rtg_dir).resolve()

    # ── 1. Load variants ─────────────────────────────────────────────
    vcf_path = rtg_dir / f"{cov}x" / caller / f"{error_type.lower()}.vcf.gz"
    print(f"\n{'='*60}")
    print(f"Loading {caller} {error_type} @ {cov}x from {vcf_path} …")
    metadata = load_vcf_rows(vcf_path, caller, error_type)
    print(f"  Variants loaded: {len(metadata)}")

    if metadata.empty:
        print("  No variants found — nothing to score.")
        return

    # ── 2. Score variants in batches ─────────────────────────────────
    unique = metadata.drop_duplicates("variant_id").reset_index(drop=True)
    n_unique = len(unique)
    batch_size = args.batch_size
    n_batches = (n_unique + batch_size - 1) // batch_size
    print(f"  Unique variants to score: {n_unique}")
    print(f"  Processing in {n_batches} batch(es) of up to {batch_size}")

    out_csv = out_dir / f"{caller}_{cov}x_{error_type}_variant_scores_tidy.csv"
    header_written = False

    pbar = tqdm(total=n_unique, desc=f"{caller} {cov}x {error_type}")
    for batch_idx in range(n_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, n_unique)
        batch = unique.iloc[start:end]

        batch_results = []
        for _, row in batch.iterrows():
            variant = genome.Variant(
                chromosome=str(row.CHROM),
                position=int(row.POS),
                reference_bases=row.REF,
                alternate_bases=row.ALT,
                name=row.variant_id,
            )
            interval = variant.reference_interval.resize(sequence_length)

            variant_scores = dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=selected_scorers,
                organism=organism,
            )
            batch_results.append(variant_scores)
            pbar.update(1)

        # ── Tidy this batch ──────────────────────────────────────────
        df_batch = variant_scorers.tidy_scores(batch_results)

        # Drop columns with unhashable objects (e.g. Variant)
        for col in list(df_batch.columns):
            if df_batch[col].dtype == "object" and col != "variant_id":
                try:
                    hash(df_batch[col].iloc[0])
                except TypeError:
                    df_batch.drop(columns=[col], inplace=True)

        # Add caller / error_type / coverage columns
        df_batch["caller"] = caller
        df_batch["error_type"] = error_type
        df_batch["coverage"] = f"{cov}x"

        # ── Append to CSV ────────────────────────────────────────────
        df_batch.to_csv(
            out_csv, mode="a", index=False, header=(not header_written)
        )
        header_written = True

        # ── Free memory ──────────────────────────────────────────────
        del batch_results, df_batch
        gc.collect()

    pbar.close()
    print(f"  Scores saved → {out_csv}")

    # ── Optional parquet ─────────────────────────────────────────────
    if not args.no_parquet:
        df = pd.read_csv(out_csv)
        pq = out_csv.with_suffix(".parquet")
        df.to_parquet(pq, index=False)
        print(f"  Parquet     → {pq}")
        del df
        gc.collect()


# ── main ─────────────────────────────────────────────────────────────
def main() -> None:
    args = parse_args()

    caller = args.caller.strip()
    error_type = args.error_type.strip().upper()

    root = args.root.resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Determine coverages ──────────────────────────────────────────
    cov_arg = args.coverage.strip().lower()
    if cov_arg == "all":
        coverages = ALL_COVERAGES
    else:
        coverages = [int(cov_arg)]

    # ── Prepare AlphaGenome model (once) ─────────────────────────────
    organism_map = {
        "human": dna_client.Organism.HOMO_SAPIENS,
        "mouse": dna_client.Organism.MUS_MUSCULUS,
    }
    organism = organism_map[args.organism]

    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f"SEQUENCE_LENGTH_{args.sequence_length}"
    ]

    api_key = args.api_key or os.getenv(args.api_key_env, "")
    assert api_key, (
        f"No API key. Use --api-key or export {args.api_key_env}=<key>"
    )
    dna_model = dna_client.create(api_key)

    # ── Select scorers ───────────────────────────────────────────────
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS

    scorer_arg = args.scorers.strip().lower()
    if scorer_arg == "all":
        selected_scorers = list(all_scorers.values())
    else:
        keys = [k.strip() for k in scorer_arg.split(",")]
        selected_scorers = [all_scorers[k] for k in keys]

    # drop unsupported for chosen organism
    unsupported = [
        s for s in selected_scorers
        if (organism.value
            not in variant_scorers.SUPPORTED_ORGANISMS[s.base_variant_scorer])
        or (s.requested_output == dna_client.OutputType.PROCAP
            and organism == dna_client.Organism.MUS_MUSCULUS)
    ]
    for s in unsupported:
        print(f"  Excluding {s} (not supported for {organism})")
        selected_scorers.remove(s)

    print(f"  Using {len(selected_scorers)} scorer(s), "
          f"organism={args.organism}, seq_len={args.sequence_length}")
    print(f"  Caller: {caller}  |  Error type: {error_type}")
    print(f"  Coverages: {coverages}")

    # ── Run for each coverage ────────────────────────────────────────
    for cov in coverages:
        score_one(args, caller, error_type, cov,
                  dna_model, organism, sequence_length,
                  selected_scorers, out_dir)

    print("\nDone ✓")


if __name__ == "__main__":
    main()
