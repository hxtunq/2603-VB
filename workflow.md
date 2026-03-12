## Workflow

Complete pipeline for multi-caller variant calling benchmarking on simulated WGS (hg38 chr22). Covers: data preparation, 5 variant callers, benchmarking with RTG vcfeval and hap.py, concordance analysis, error pattern extraction, stratified annotation, and visualization.

Each runnable shell block starts with a `# pwd:` comment showing the expected working directory.

### Repository Structure

```
config/config.sh                          # Central configuration
pipelines/
  _sentieon_common.sh                     # Sentieon utility functions
  03_call_hc.sh                           # GATK HaplotypeCaller
  04_call_dv.sh                           # DeepVariant (Docker)
  05_call_strelka.sh                      # Strelka2 (Docker)
  06_call_freebayes.sh                    # FreeBayes
  07_call_dnascope.sh                     # Sentieon DNAscope (shared BAM)
  07_call_dnascope_fastq.sh              # Sentieon DNAscope (raw FASTQ)
evaluation/
  _happy_common.sh                        # hap.py utility functions
  eval_happy.sh                           # hap.py evaluation (Docker)
  eval_vcfeval.sh                         # RTG vcfeval evaluation
  gather_stats.sh                         # Aggregate hap.py CSVs → TSV
  build_variant_table.sh                  # Binary concordance matrix from VCFs
  generate_gc_strata.py                   # GC% BED strata from reference
  generate_cds_gc_strata.py               # GC% BED strata within CDS
  prepare_stratification_hg38.sh          # CDS + coverage + GIAB strata
  prepare_cds_gc_stratification.sh        # CDS-specific GC strata wrapper
  prepare_clinvar.sh                      # ClinVar pathogenic BED
  prepare_giab_strata.sh                  # GIAB repeat/mappability/segdup BEDs
  concordance/
    run_concordance.sh                    # Pairwise RTG vcfeval
    concordance_matrix.py                 # Parse pairwise → concordance TSV
downstream/
  build_concordance_from_vcfs.py          # Rebuild concordance tables from VCFs
  extract_error_patterns.py               # Extract FP/FN error patterns
  annotate_variants_strata.py             # Annotate variants with genomic context
visualization/
  benchmark_plots.R                       # Per-coverage F1/Precision/Recall bars
  plot_summary.py                         # Cross-coverage heatmaps, grouped bars
  plot_coverage_comparison.py             # Coverage-impact line/radar/area plots
  plot_concordance.py                     # Caller-vs-caller concordance heatmaps
  01_pairwise_wilcoxon.R                  # Pairwise Wilcoxon signed-rank heatmaps
  02_pca_scatterplot.R                    # PCA on discordant variants
  03_stratified_plots.R                   # Stratified F1 + genomic context heatmap
  04_clinvar_analysis.R                   # ClinVar detection rate analysis
  upset-plot.R                            # UpSet plots for FP/FN overlap
```

## 1. Setup

```bash
git clone <repo-url> variant-calling-benchmark
cd variant-calling-benchmark

mkdir -p data/reference
mkdir -p data/simulated
mkdir -p results/preprocessing
mkdir -p results/variants
mkdir -p results/eval
mkdir -p results/plots
mkdir -p logs

source config/config.sh
```

## 2. Reference Assets

### 2.1 Reference FASTA

```bash
# pwd: variant-calling-benchmark/data/reference

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz

samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

### 2.2 Known Sites For BQSR

```bash
# pwd: variant-calling-benchmark/data/reference

# dbSNP
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.vcf.gz
bcftools view -r chr22 -Oz -o dbsnp138.hg38.chr22.vcf.gz dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.chr22.vcf.gz

# Mills & 1000G gold-standard indels
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
bcftools view -r chr22 -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

# 1000G high-confidence SNPs
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
bcftools view -r chr22 -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz \
  1000G_phase1.snps.high_confidence.hg38.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
```

### 2.3 BED Assets And RTG SDF

```bash
# pwd: variant-calling-benchmark/data/reference

# High-confidence BED: non-N intervals ≥1000bp (requires Biopython)
python3 - <<'PY'
from Bio import SeqIO
import re

with open("chr22_highconf.bed", "w") as out:
    for record in SeqIO.parse("chr22.fa", "fasta"):
        for match in re.finditer(r"[ACGT]+", str(record.seq).upper()):
            if match.end() - match.start() >= 1000:
                out.write(f"{record.id}\t{match.start()}\t{match.end()}\n")
PY

# RTG SDF template (required for vcfeval and hap.py --engine vcfeval)
rtg format -o chr22.sdf chr22.fa
```

> **Note:** `stratification_chr22.tsv` is generated later by `evaluation/prepare_stratification_hg38.sh`.

### 2.4 Sentieon Assets (Optional — DNAscope Only)

```bash
# pwd: variant-calling-benchmark

# 1. License & software
mkdir -p .tools/sentieon
cp /path/to/Hanoi_University_Of_Science_And_Technology_eval.lic .tools/sentieon/

wget -O .tools/sentieon/sentieon-genomics-202503.02.tar.gz \
  "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-202503.02.tar.gz"
tar -xzf .tools/sentieon/sentieon-genomics-202503.02.tar.gz -C .tools/sentieon/

# 2. Install sentieon-cli
pip install .tools/sentieon-cli-src/

# 3. Environment variables
export SENTIEON_LICENSE="$PWD/.tools/sentieon/Hanoi_University_Of_Science_And_Technology_eval.lic"
export SENTIEON_BIN_DIR="$PWD/.tools/sentieon/sentieon-genomics-202503.02/bin"
export PATH="${SENTIEON_BIN_DIR}:$PATH"

# Verify
command -v sentieon-cli
command -v sentieon

# 4. DNAscope model bundle
mkdir -p data/reference/models
wget -O data/reference/models/SentieonIlluminaWGS2.2.bundle \
  https://s3.amazonaws.com/sentieon-release/other/SentieonIlluminaWGS2.2.bundle
```

### 2.5 Stratification BEDs

All stratification is handled by dedicated scripts. Run them in this order:

```bash
# pwd: variant-calling-benchmark

# GC content strata (7 bins: 0–20%, 20–30%, ..., 80–100%)
python evaluation/generate_gc_strata.py data/reference/chr22.fa data/reference/gc_strata 1000

# GENCODE CDS + CDS vicinity + coverage strata + GIAB repeat/mappability/segdup
bash evaluation/prepare_stratification_hg38.sh

# %GC within CDS regions (requires CDS-canonical BED from above)
bash evaluation/prepare_cds_gc_stratification.sh

# ClinVar pathogenic variants on chr22
bash evaluation/prepare_clinvar.sh
```

**`prepare_stratification_hg38.sh`** performs 4 steps:
1. Downloads GENCODE v44 GTF → extracts protein-coding CDS and canonical CDS BEDs for chr22
2. Generates CDS vicinity BEDs (±0, ±25, ±50, ±100 bp padding)
3. Computes coverage-stratification BEDs (normalized depth bins per BAM)
4. Calls `prepare_giab_strata.sh` → downloads GIAB v3.3 BEDs:
   - Tandem repeats + homopolymers (`GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed`)
   - Low mappability (`GRCh38_lowmappabilityall.bed`)
   - Segmental duplications >10kb (`GRCh38_segdups_gt10kb.bed`)

All BED paths are assembled into `data/reference/stratification_chr22.tsv` for hap.py.

## 3. Truth Set And Simulated FASTQs

### 3.1 Truth VCF From simuG

```bash
# pwd: variant-calling-benchmark/data

git clone https://github.com/yjx1217/simuG.git

perl simuG/simuG.pl \
  -refseq reference/chr22.fa \
  -snp_count 18000 \
  -indel_count 2500 \
  -titv_ratio 2.0 \
  -prefix simulated/SIMULATED_SAMPLE_chr22

# Reheader and merge SNP + INDEL VCFs
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.vcf
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.vcf

bcftools concat \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf | \
bcftools sort -Oz -o simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

tabix -p vcf simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz
```

### 3.2 WGS FASTQs (ART Illumina)

```bash
# pwd: variant-calling-benchmark/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

for COV in 10 20 30 50; do
  art_illumina \
    -ss HS25 -i "${MUTATED_FASTA}" -p -l 150 -f "${COV}" \
    -m 350 -s 50 -rs 42 -o "${SIM_DIR}/${PREFIX}_${COV}x_" -na

  mv "${SIM_DIR}/${PREFIX}_${COV}x_1.fq" "${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq"
  mv "${SIM_DIR}/${PREFIX}_${COV}x_2.fq" "${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq"
done
```

## 4. Shared Preprocessing

### 4.1 Alignment And Sort

```bash
# pwd: variant-calling-benchmark
PREFIX="SIMULATED_SAMPLE_chr22"
REF="data/reference/chr22.fa"

for COV in 10 20 30 50; do
  R1="data/simulated/${PREFIX}_${COV}x_R1.fastq.gz"
  R2="data/simulated/${PREFIX}_${COV}x_R2.fastq.gz"
  OUTDIR="results/preprocessing/${COV}x"
  mkdir -p "${OUTDIR}" "logs/${COV}x"

  bwa mem -t 4 -M \
    -R "@RG\tID:${PREFIX}_${COV}x\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${R1}" "${R2}" | \
    samtools sort -@ 4 -m 2G -o "${OUTDIR}/${PREFIX}_aligned.bam" -

  samtools index "${OUTDIR}/${PREFIX}_aligned.bam"
done
```

### 4.2 MarkDuplicates And Coverage Stats

```bash
# pwd: variant-calling-benchmark
PREFIX="SIMULATED_SAMPLE_chr22"

for COV in 10 20 30 50; do
  OUTDIR="results/preprocessing/${COV}x"

  gatk MarkDuplicates \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_dedup.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true

  samtools stats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
  samtools flagstat "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
  samtools idxstats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_idxstats.txt"
  mosdepth -t 4 --by 1000 "${OUTDIR}/${PREFIX}_coverage" "${OUTDIR}/${PREFIX}_dedup.bam"
done
```

## 5. Variant Calling (5 Callers × 4 Coverages)

```bash
# pwd: variant-calling-benchmark
for COV in 10 20 30 50; do
  bash pipelines/03_call_hc.sh "${COV}"           # GATK HaplotypeCaller (BQSR + hard filter)
  bash pipelines/04_call_dv.sh "${COV}"            # DeepVariant (Docker)
  bash pipelines/05_call_strelka.sh "${COV}"       # Strelka2 (Docker)
  bash pipelines/06_call_freebayes.sh "${COV}"     # FreeBayes (+ GATK left-align/filter)
  bash pipelines/07_call_dnascope.sh "${COV}"      # Sentieon DNAscope (shared BAM)
done
```

Optional: DNAscope from raw FASTQ (includes alignment internally):

```bash
for COV in 10 20 30 50; do
  bash pipelines/07_call_dnascope_fastq.sh "${COV}"
done
```

**Output per caller/coverage:** `results/variants/{COV}x/{caller}/SS_chr22_{TAG}_{COV}x_WGS.vcf.gz`

| Script | Caller | Tag | Notes |
|---|---|---|---|
| `03_call_hc.sh` | GATK HaplotypeCaller | `HC` | BQSR → call → hard-filter SNP+INDEL |
| `04_call_dv.sh` | DeepVariant 1.9.0 | `DV` | Docker; also produces gVCF |
| `05_call_strelka.sh` | Strelka2 2.9.10 | `STRELKA` | Docker; configure + run in 2 steps |
| `06_call_freebayes.sh` | FreeBayes 1.3.10 | `FB` | Raw call → left-align → split → filter |
| `07_call_dnascope.sh` | DNAscope (shared BAM) | `DNASCOPE` | Requires Sentieon license |
| `07_call_dnascope_fastq.sh` | DNAscope (FASTQ) | `DNASCOPE_FASTQ` | FASTQ → align → call (one CLI) |

## 6. Benchmarking

### 6.1 RTG vcfeval

```bash
# pwd: variant-calling-benchmark
bash evaluation/eval_vcfeval.sh
```

Evaluates all 5 callers × 4 coverages against truth VCF using RTG vcfeval.

**Outputs:**
- Per-caller: `results/eval/vcfeval/{COV}x/{caller}/` → `tp.vcf.gz`, `fp.vcf.gz`, `fn.vcf.gz`, `summary.txt`
- Aggregated: `results/eval/vcfeval_all_stats.tsv`

### 6.2 hap.py (Docker)

```bash
# pwd: variant-calling-benchmark
bash evaluation/eval_happy.sh
```

Evaluates using hap.py with stratification (requires Docker). Calls `gather_stats.sh` automatically.

**Outputs:**
- Per-caller: `results/eval/{COV}x/WGS/{caller}_eval_data/report.extended.csv`
- Aggregated: `results/eval/all_stats.tsv` (with columns: Coverage, Mode, Caller, Type, Filter, Subtype, Subset, METRIC.\*, TRUTH.\*, QUERY.\*)

### 6.3 Concordance Analysis

```bash
# pwd: variant-calling-benchmark

# Pairwise RTG vcfeval (all 10 caller pairs × 4 coverages)
bash evaluation/concordance/run_concordance.sh

# Parse pairwise results → concordance matrix TSV
python evaluation/concordance/concordance_matrix.py \
    results/eval/concordance \
    results/eval/concordance/concordance_matrix.tsv
```

**Output:** `concordance_matrix.tsv` — columns: `Coverage, CallerA, CallerB, Concordance`

### 6.4 Variant-Level Concordance Tables

```bash
# pwd: variant-calling-benchmark

# Build binary concordance matrix from vcfeval TP/FP/FN VCFs
bash evaluation/build_variant_table.sh
```

**Outputs per coverage:**
- `results/eval/vcfeval/variant_concordance_{COV}x.tsv` — binary matrix (variant_id, truth, HC, DV, ST, FB, DS)
- `results/eval/vcfeval/upset_long_callset_{COV}x.csv` — FP long format for UpSet
- `results/eval/vcfeval/upset_long_baseline_{COV}x.csv` — FN long format for UpSet

Alternative (Python, if TSV files are missing):

```bash
python downstream/build_concordance_from_vcfs.py \
    --project-root . --search-root results/eval \
    --out-dir results/eval/vcfeval --coverages 10,20,30,50
```

## 7. Downstream Error Analysis

### 7.1 Extract Error Patterns

```bash
# pwd: variant-calling-benchmark
python downstream/extract_error_patterns.py \
    --project-root . \
    --vcfeval-dir results/eval/vcfeval \
    --out-dir results/analysis/error_patterns \
    --coverages 10,20,30,50 \
    --callers HC,DV,ST,FB,DS
```

Extracts 6 error pattern types per coverage:

| Pattern | Description |
|---|---|
| `all_5_fn` | All 5 callers missed (truth=1, n_called=0) |
| `all_5_fp` | All 5 callers false positive (truth=0, n_called=5) |
| `single_caller_fn` | Exactly 1 caller missed, other 4 detected |
| `single_caller_fp` | Exactly 1 caller FP, others did not call |
| `dv_ds_fn` | Both DeepVariant and DNAscope missed |
| `dv_ds_fp` | Both DeepVariant and DNAscope false positive |

**Outputs:**
- `results/analysis/error_patterns/{pattern}_{COV}x.tsv` (24 files)
- `results/analysis/error_patterns/pattern_counts_summary.tsv`

### 7.2 Annotate Variants With Genomic Context

```bash
# pwd: variant-calling-benchmark
python downstream/annotate_variants_strata.py \
    --project-root . \
    --patterns-dir results/analysis/error_patterns \
    --gc-dir data/reference/gc_strata \
    --cds-bed data/reference/stratification/CDS-canonical.chr22.bed \
    --repeats-dir data/reference/repeats \
    --out-dir results/analysis/stratification
```

Annotates each error-pattern variant with genomic region overlaps:

| Column | Source BED | Description |
|---|---|---|
| `gc_bin` | GC strata | GC content bin (e.g. GC_30_40) |
| `gc_pct_mid` | (derived) | Midpoint of GC bin (e.g. 35.0) |
| `cds_gc_bin` | CDS GC strata | GC bin within CDS regions |
| `in_cds_canonical` | CDS-canonical BED | Overlaps Ensembl canonical CDS (0/1) |
| `in_protein_coding` | CDS BED | Overlaps protein-coding CDS (0/1) |
| `in_tandem_repeat` | GIAB v3.3 | Overlaps tandem repeat or homopolymer (0/1) |
| `in_low_mapq` | GIAB v3.3 | Overlaps low-mappability region (0/1) |
| `in_segdup` | GIAB v3.3 | Overlaps segmental duplication >10kb (0/1) |

**Outputs:**
- `results/analysis/stratification/{pattern}_{COV}x_with_strata.tsv` (24 files)
- `results/analysis/stratification/strata_counts_summary.tsv`

## 8. Visualization

### 8.1 Summary Plots

```bash
# pwd: variant-calling-benchmark

# R: per-coverage F1/Precision/Recall bar charts
Rscript visualization/benchmark_plots.R results/eval/all_stats.tsv results/plots

# Python: cross-coverage heatmaps, grouped bars, SNP vs INDEL breakdown
python visualization/plot_summary.py results/eval/all_stats.tsv results/plots
```

**R outputs per (Mode, Coverage):** `fig_f1_*.png`, `fig_precision_recall_*.png`, `fig_counts_*.png`

**Python cross-coverage outputs:** `grouped_bars_snp.png`, `grouped_bars_indel.png`, `f1_heatmap_wgs.png`, `snp_vs_indel_wgs.png`, `variant_breakdown_wgs.png`

### 8.2 Coverage Comparison Plots

```bash
# pwd: variant-calling-benchmark
python visualization/plot_coverage_comparison.py \
    results/eval/vcfeval_all_stats.tsv results/plots
```

**Outputs:** `coverage_lines.png`, `coverage_delta_bars.png`, `coverage_radar.png`, `coverage_fn_fp_area.png`, `coverage_f1_heatmap.png`

### 8.3 Concordance Heatmaps

```bash
# pwd: variant-calling-benchmark
python visualization/plot_concordance.py \
    results/eval/concordance/concordance_matrix.tsv results/plots
```

**Outputs:** `concordance_heatmap_{COV}.png` (per coverage), `concordance_heatmap_average.png`

### 8.4 Advanced Analysis

```bash
# pwd: variant-calling-benchmark

# Pairwise Wilcoxon signed-rank test heatmaps (statistical significance)
Rscript visualization/01_pairwise_wilcoxon.R results/eval/all_stats.tsv results/plots

# PCA on discordant variants (caller clustering)
Rscript visualization/02_pca_scatterplot.R

# Stratified F1 + genomic context heatmap/bar chart
Rscript visualization/03_stratified_plots.R results/eval/all_stats.tsv results/plots

# ClinVar pathogenic variant detection analysis
Rscript visualization/04_clinvar_analysis.R results/eval/all_stats.tsv results/plots

# UpSet plots (FP/FN overlap across callers)
Rscript visualization/upset-plot.R
```

**`03_stratified_plots.R`** generates 8 plot types:

| # | Plot | Description |
|---|---|---|
| 1 | `boxplot_f1_by_coverage_*.png` | F1 distribution per caller × coverage |
| 2 | `lineplot_f1_by_gc_*.png` | F1 vs %GC content per caller |
| 3 | `lineplot_f1_by_covstratum_*.png` | F1 vs normalized coverage depth |
| 4 | `lineplot_f1_by_cds_distance_*.png` | F1 vs distance from CDS boundary |
| 5 | `boxplot_f1_cds_comparison_*.png` | F1 in Overall vs CDS vs Canonical CDS |
| 6 | `lineplot_f1_by_cds_gc_*.png` | F1 vs %GC within CDS regions |
| 7 | `heatmap_caller_vs_region.png` | Caller × Region → % FP/FN errors (from downstream strata) |
| 8 | `bar_errors_by_region.png` | Error counts per difficult region, grouped by caller |

Plots 7–8 require `results/analysis/stratification/*_with_strata.tsv` (produced by step 7.2).

## 9. Complete Pipeline — Quick Reference

```bash
# pwd: variant-calling-benchmark
source config/config.sh

# ── Reference & stratification ──
python evaluation/generate_gc_strata.py data/reference/chr22.fa data/reference/gc_strata 1000
bash evaluation/prepare_stratification_hg38.sh
bash evaluation/prepare_cds_gc_stratification.sh
bash evaluation/prepare_clinvar.sh

# ── Calling ──
for COV in 10 20 30 50; do
  bash pipelines/03_call_hc.sh "${COV}"
  bash pipelines/04_call_dv.sh "${COV}"
  bash pipelines/05_call_strelka.sh "${COV}"
  bash pipelines/06_call_freebayes.sh "${COV}"
  bash pipelines/07_call_dnascope.sh "${COV}"
done

# ── Benchmarking ──
bash evaluation/eval_vcfeval.sh
bash evaluation/eval_happy.sh
bash evaluation/concordance/run_concordance.sh
python evaluation/concordance/concordance_matrix.py \
    results/eval/concordance results/eval/concordance/concordance_matrix.tsv
bash evaluation/build_variant_table.sh

# ── Downstream analysis ──
python downstream/extract_error_patterns.py
python downstream/annotate_variants_strata.py

# ── Visualization ──
Rscript visualization/benchmark_plots.R results/eval/all_stats.tsv results/plots
python visualization/plot_summary.py results/eval/all_stats.tsv results/plots
python visualization/plot_coverage_comparison.py results/eval/vcfeval_all_stats.tsv results/plots
python visualization/plot_concordance.py results/eval/concordance/concordance_matrix.tsv results/plots
Rscript visualization/01_pairwise_wilcoxon.R
Rscript visualization/02_pca_scatterplot.R
Rscript visualization/03_stratified_plots.R
Rscript visualization/04_clinvar_analysis.R
Rscript visualization/upset-plot.R
```

## 10. Method Notes

- **Coverages:** 10x / 20x / 30x / 50x — spans low-depth research to production-grade WGS.
- **Callers:** GATK HC, DeepVariant, Strelka2, FreeBayes, Sentieon DNAscope — represents rule-based, ML-based, and hybrid approaches.
- **Evaluation:** Dual-engine (RTG vcfeval + hap.py) ensures robust benchmarking with stratification support.
- **Stratification layers:** GC content (7 bins), CDS regions (canonical + vicinity), coverage depth (9 bins), CDS-specific GC, tandem repeats & homopolymers (GIAB v3.3), low mappability (GIAB v3.3), segmental duplications (GIAB v3.3), ClinVar pathogenic variants.
- **Error analysis:** 6 error patterns identify systematic caller weaknesses (universal FN/FP, single-caller errors, DV-DS shared errors), annotated with genomic context for region-specific failure diagnosis.

## 11. Dependencies

| Tool | Version | Purpose |
|---|---|---|
| BWA | ≥0.7.17 | Read alignment |
| Samtools | ≥1.17 | BAM handling, depth, indexing |
| GATK | 4.6.2.0 | HaplotypeCaller, BQSR, MarkDuplicates, filtering |
| FreeBayes | 1.3.10 | Bayesian variant calling |
| DeepVariant | 1.9.0 | DL-based variant calling (Docker) |
| Strelka2 | 2.9.10 | Germline variant calling (Docker) |
| Sentieon | 202503.02 | DNAscope variant calling (license required) |
| RTG Tools | ≥3.12.1 | vcfeval benchmarking |
| hap.py | v0.3.12 | GA4GH benchmarking (Docker) |
| bedtools | ≥2.30 | BED interval operations |
| bcftools | ≥1.17 | VCF manipulation |
| mosdepth | ≥0.3 | Coverage statistics |
| ART | — | Illumina read simulation |
| simuG | — | Variant simulation |
| Python | ≥3.8 | pandas, downstream scripts |
| R | ≥4.0 | ggplot2, dplyr, tidyr, cowplot, stringr, ComplexHeatmap, UpSetR |