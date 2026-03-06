# Variant Calling Benchmark

Benchmark variant callers on simulated `chr22` Illumina data with one primary comparison track and one optional Sentieon end-to-end path:

- Primary benchmark: fair head-to-head comparison on the same shared dedup BAM
- Optional DNAscope FASTQ runs: end-to-end Sentieon alignment + calling kept separate from the fair comparison

The repository keeps WGS and WES as separate study modes.

## Current Scope

Implemented benchmark components:

- Shared preprocessing: alignment, sorting, duplicate marking, coverage stats
- Shared-BAM callers:
  - GATK HaplotypeCaller
  - DeepVariant
  - Strelka2 + Manta
  - FreeBayes
  - Sentieon DNAscope
- Optional end-to-end DNAscope from raw FASTQ:
  - WGS
  - WES
- Truth-based evaluation with hap.py + RTG vcfeval
- Summary plots from `results/eval/all_stats.tsv`

Not implemented in this repo:

- SnpEff / dbNSFP / ClinVar annotation pipeline
- AlphaGenome functional-risk scoring

Those items are future work only, not runnable workflow steps.

## Study Design

Main coverage grid:

- WGS: `10x 20x 30x 50x`
- WES: `50x 100x 200x`

Design rules:

- The shared-BAM track is the primary benchmark result because every caller sees the same dedup BAM.
- The optional DNAscope FASTQ runs are supplementary because alignment is part of the pipeline and cannot be interpreted as a fair caller-only comparison.
- The archived `-Impact-of-Sequencing-Depth-on-Variant-Detection--main` project is an HG002 chr22 exome downsampling study. Its findings are useful only as qualitative coverage guidance, not as WGS-vs-WES evidence.

## Repository Layout

```text
variant-calling-benchmark/
├── config/
│   └── config.sh
├── evaluation/
│   ├── _happy_common.sh
│   ├── eval_happy.sh
│   └── gather_stats.sh
├── pipelines/
│   ├── _sentieon_common.sh
│   ├── 03_call_hc*.sh
│   ├── 04_call_dv*.sh
│   ├── 05_call_strelka*.sh
│   ├── 06_call_freebayes*.sh
│   ├── 07_call_dnascope*.sh
│   └── 07_call_dnascope_fastq*.sh
├── visualization/
│   ├── benchmark_plots.R
│   └── plot_summary.py
├── workflow.md
└── archive/
```

## Quick Start

`workflow.md` contains the detailed setup commands. The high-level execution order is:

```bash
source config/config.sh

# Build RTG SDF once before hap.py/vcfeval
rtg format -o "${RTG_SDF}" "${REF_FASTA}"

# Prepare references, truth set, simulated FASTQs, and shared preprocessing
# See workflow.md for the full commands

# Shared-BAM benchmark
for COV in 10 20 30 50; do
  bash pipelines/03_call_hc.sh "${COV}"
  bash pipelines/04_call_dv.sh "${COV}"
  bash pipelines/05_call_strelka.sh "${COV}"
  bash pipelines/06_call_freebayes.sh "${COV}"
  bash pipelines/07_call_dnascope.sh "${COV}"
done

for COV in 50 100 200; do
  bash pipelines/03_call_hc_wes.sh "${COV}"
  bash pipelines/04_call_dv_wes.sh "${COV}"
  bash pipelines/05_call_strelka_wes.sh "${COV}"
  bash pipelines/06_call_freebayes_wes.sh "${COV}"
  bash pipelines/07_call_dnascope_wes.sh "${COV}"
done

# Optional: end-to-end DNAscope from raw FASTQs
for COV in 10 20 30 50; do
  bash pipelines/07_call_dnascope_fastq.sh "${COV}"
done

for COV in 50 100 200; do
  bash pipelines/07_call_dnascope_fastq_wes.sh "${COV}"
done

# Evaluation
bash evaluation/eval_happy.sh

# Visualization
Rscript visualization/benchmark_plots.R results/eval/all_stats.tsv results/plots
python visualization/plot_summary.py results/eval/all_stats.tsv results/plots
```

## Outputs

Primary benchmark:

- Variant calls: `results/variants/{coverage}/...` and `results/variants/{coverage}_wes/...`
- Evaluation: `results/eval/...`
- Aggregated stats: `results/eval/all_stats.tsv`

Benchmark runtime / CPU / RSS are appended per coverage to:

- `logs/{coverage}/benchmark_metrics.tsv`
- `logs/{coverage}_wes/benchmark_metrics.tsv`

## Sentieon Notes

Shared-BAM DNAscope and raw-FASTQ DNAscope use the local `sentieon-cli` command. Keep `sentieon-cli` in `PATH`.

Required inputs for DNAscope:

- `SENTIEON_LICENSE`
- `DNASCOPE_WGS_MODEL`
- `DNASCOPE_WES_MODEL`

The shared-BAM DNAscope scripts use the shared dedup BAM directly. This is the fair-comparison path and matches official Sentieon support for variant calling from sorted BAM/CRAM.

The optional `07_call_dnascope_fastq*.sh` scripts run the full DNAscope alignment + calling pipeline from raw FASTQ for WGS and WES. Treat those as separate end-to-end experiments, not as shared-BAM fair-comparison runs.

## References

- Sentieon CLI docs: https://support.sentieon.com/docs/sentieon_cli/
- Sentieon models: https://github.com/Sentieon/sentieon-models
- Archived depth study: `archive/-Impact-of-Sequencing-Depth-on-Variant-Detection--main/`
