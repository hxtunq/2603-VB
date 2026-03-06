# Variant Calling Benchmark

Benchmark variant callers on simulated `chr22` Illumina data with two explicit tracks:

- Track A: fair head-to-head comparison on the same shared dedup BAM
- Track B: supplementary Sentieon pangenome end-to-end pipeline, reported separately

The repository keeps WGS and WES as separate study modes and does not merge pangenome results into the main ranking.

## Current Scope

Implemented benchmark components:

- Shared preprocessing: alignment, sorting, duplicate marking, coverage stats
- Track A callers:
  - GATK HaplotypeCaller
  - DeepVariant
  - Strelka2 + Manta
  - FreeBayes
  - Sentieon DNAscope
- Track B callers:
  - Sentieon DNAscope Pangenome WGS
  - Sentieon DNAscope Pangenome WES
- Truth-based evaluation with hap.py + RTG vcfeval
- Summary plots from `results/eval/all_stats.tsv` and `results/eval_track_b/all_stats.tsv`

Not implemented in this repo:

- SnpEff / dbNSFP / ClinVar annotation pipeline
- AlphaGenome functional-risk scoring

Those items are now future work only, not runnable workflow steps.

## Study Design

Main coverage grid:

- WGS: `10x 20x 30x 50x`
- WES: `50x 100x 200x`

Design rules:

- Track A is the primary benchmark result because every caller sees the same dedup BAM.
- Track B is supplementary because pangenome alignment changes the input BAM and cannot be interpreted as a fair caller-only comparison.
- The archived `-Impact-of-Sequencing-Depth-on-Variant-Detection--main` project is an HG002 chr22 exome downsampling study. Its findings are useful only as qualitative coverage guidance, not as WGS-vs-WES evidence.

## Repository Layout

```text
variant-calling-benchmark/
├── config/
│   └── config.sh
├── evaluation/
│   ├── _happy_common.sh
│   ├── eval_happy.sh
│   ├── eval_happy_track_b.sh
│   └── gather_stats.sh
├── pipelines/
│   ├── _sentieon_common.sh
│   ├── 03_call_hc*.sh
│   ├── 04_call_dv*.sh
│   ├── 05_call_strelka*.sh
│   ├── 06_call_freebayes*.sh
│   ├── 07_call_dnascope*.sh
│   └── 07_call_dnascope_pangenome*.sh
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

# Track A: fair comparison on shared BAMs
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

# Track B: supplementary Sentieon pangenome runs
for COV in 10 20 30 50; do
  bash pipelines/07_call_dnascope_pangenome.sh "${COV}"
done

for COV in 50 100 200; do
  bash pipelines/07_call_dnascope_pangenome_wes.sh "${COV}"
done

# Evaluation
bash evaluation/eval_happy.sh
bash evaluation/eval_happy_track_b.sh

# Visualization
Rscript visualization/benchmark_plots.R results/eval/all_stats.tsv results/plots
python visualization/plot_summary.py results/eval/all_stats.tsv results/plots
```

## Outputs

Track A:

- Variant calls: `results/variants/{coverage}/...` and `results/variants/{coverage}_wes/...`
- Evaluation: `results/eval/...`
- Aggregated stats: `results/eval/all_stats.tsv`

Track B:

- Variant calls: `results/variants/{coverage}/dnascope_pangenome/...` and `results/variants/{coverage}_wes/dnascope_pangenome_wes/...`
- Evaluation: `results/eval_track_b/...`
- Aggregated stats: `results/eval_track_b/all_stats.tsv`

Benchmark runtime / CPU / RSS are appended per coverage to:

- `logs/{coverage}/benchmark_metrics.tsv`
- `logs/{coverage}_wes/benchmark_metrics.tsv`

## Sentieon Notes

The repo is standardized on Docker for Sentieon execution.

Required inputs for DNAscope:

- `SENTIEON_LICENSE`
- `DNASCOPE_WGS_MODEL`
- `DNASCOPE_WES_MODEL`
- `PANGENOME_INDEX` for Track B

Track A DNAscope uses the shared dedup BAM directly. This is the fair-comparison path and matches official Sentieon support for variant calling from sorted BAM/CRAM.

Track B DNAscope Pangenome starts from the simulated FASTQs and uses `bwa-mem2-pangenome`, so its result must stay separate from the main benchmark tables.

## References

- Sentieon pangenome usage: https://support.sentieon.com/versions/202503.02/docs/Pangenome_usage/pangenome/
- Sentieon CLI docs: https://support.sentieon.com/docs/sentieon_cli/
- Sentieon models: https://github.com/Sentieon/sentieon-models
- Archived depth study: `archive/-Impact-of-Sequencing-Depth-on-Variant-Detection--main/`
