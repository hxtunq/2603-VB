# Variant Calling Benchmark

Benchmark variant callers on simulated `chr22` Illumina WGS data with one primary comparison track and one optional Sentieon end-to-end path:

- Primary benchmark: fair head-to-head comparison on the same shared dedup BAM
- Optional DNAscope FASTQ runs: end-to-end Sentieon alignment + calling kept separate from the fair comparison

## Current Scope

Implemented benchmark components:

- Shared preprocessing: alignment, sorting, duplicate marking, coverage stats
- Shared-BAM callers:
  - GATK HaplotypeCaller
  - DeepVariant
  - Strelka2 + Manta
  - FreeBayes
  - Sentieon DNAscope
- Optional end-to-end DNAscope from raw FASTQ
- Truth-based evaluation with hap.py + RTG vcfeval
- Summary plots from `results/eval/all_stats.tsv`
- Advanced analysis: Pairwise Wilcoxon heatmaps, PCA, Stratified performance (CDS/GC/Coverage), UpSet plots, ClinVar pathogenic detection

Not implemented in this repo:

- SnpEff / dbNSFP annotation pipeline
- AlphaGenome functional-risk scoring

Those items are future work only, not runnable workflow steps.

## Study Design

Coverage grid: `10x 20x 30x 50x`

Design rules:

- The shared-BAM track is the primary benchmark result because every caller sees the same dedup BAM.
- The optional DNAscope FASTQ runs are supplementary because alignment is part of the pipeline and cannot be interpreted as a fair caller-only comparison.

## Repository Layout

```text
variant-calling-benchmark/
├── config/
│   └── config.sh
├── evaluation/
│   ├── _happy_common.sh
│   ├── eval_happy.sh
│   ├── fix_truth_vcf.sh
│   ├── gather_stats.sh
│   ├── generate_gc_strata.py        # GC-content stratification BEDs
│   └── concordance/
│       ├── run_concordance.sh        # Pairwise RTG vcfeval
│       └── concordance_matrix.py     # Build concordance matrix
├── pipelines/
│   ├── _sentieon_common.sh
│   ├── 03_call_hc.sh
│   ├── 04_call_dv.sh
│   ├── 05_call_strelka.sh
│   ├── 06_call_freebayes.sh
│   ├── 07_call_dnascope.sh
│   └── 07_call_dnascope_fastq.sh
├── visualization/
│   ├── benchmark_plots.R
│   ├── plot_summary.py               # Per-coverage + cross-coverage plots
│   └── plot_concordance.py           # Concordance heatmaps
├── workflow.md
└── archive/
```

## Outputs

Primary benchmark:

- Variant calls: `results/variants/{coverage}/...`
- Evaluation: `results/eval/...`
- Aggregated stats: `results/eval/all_stats.tsv`
- Concordance matrix: `results/eval/concordance/concordance_matrix.tsv`
- Variant breakdown (SNP/INDEL): `results/plots/variant_breakdown_*.tsv`

Visualization outputs:

- F1 heatmap (callers × coverage): `results/plots/f1_heatmap_*.png`
- Grouped bar charts: `results/plots/grouped_bars_*.png`
- SNP vs INDEL comparison: `results/plots/snp_vs_indel_*.png`
- Concordance heatmaps: `results/plots/concordance_heatmap_*.png`

Benchmark runtime / CPU / RSS are appended per coverage to:

- `logs/{coverage}/benchmark_metrics.tsv`

## Sentieon Notes

Required inputs for DNAscope:

- `SENTIEON_LICENSE`
- `DNASCOPE_WGS_MODEL` pointing to a `.bundle` file

If you previously extracted the Sentieon model bundles and exported `DNASCOPE_WGS_MODEL` as a directory, update it to the original `.bundle` archive path instead. The scripts now resolve a sibling `.bundle` automatically when possible, but the canonical configuration is the archive file.

The shared-BAM DNAscope scripts use the shared dedup BAM directly. This is the fair-comparison path and matches official Sentieon support for variant calling from sorted BAM/CRAM.

## References

- Sentieon CLI docs: https://support.sentieon.com/docs/sentieon_cli/
- Sentieon models: https://github.com/Sentieon/sentieon-models