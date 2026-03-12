## Workflow

This file documents the implemented benchmark only. Annotation-heavy stages such as SnpEff, dbNSFP, ClinVar, and AlphaGenome are future work and are not part of the runnable pipeline in this repo.
Each runnable shell block starts with a `# pwd:` comment showing the expected working directory.

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

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.vcf.gz
bcftools view -r chr22 -Oz -o dbsnp138.hg38.chr22.vcf.gz dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.chr22.vcf.gz

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
bcftools view -r chr22 -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
bcftools view -r chr22 -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz \
  1000G_phase1.snps.high_confidence.hg38.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
```

### 2.3 BED Assets And hap.py Support Files

```bash
# pwd: variant-calling-benchmark/data/reference

# phải dùng đến env biopython: conda activate biopython

python3 - <<'PY'
from Bio import SeqIO
import re

with open("chr22_highconf.bed", "w") as out:
    for record in SeqIO.parse("chr22.fa", "fasta"):
        for match in re.finditer(r"[ACGT]+", str(record.seq).upper()):
            if match.end() - match.start() >= 1000:
                out.write(f"{record.id}\t{match.start()}\t{match.end()}\n")
PY

# quay về env ban đầu

> **Note:** The `stratification_chr22.tsv` will be generated later by `evaluation/prepare_stratification_hg38.sh`.

rtg format -o chr22.sdf chr22.fa
```

### 2.5 GC Content Stratification BEDs

```bash
# pwd: variant-calling-benchmark

python evaluation/generate_gc_strata.py \
    data/reference/chr22.fa \
    data/reference/gc_strata \
    1000
```

This generates 7 BED files (GC_0_20, GC_20_30, ..., GC_80_100) in `data/reference/gc_strata/` and updates `stratification_chr22.tsv` with the new GC bins. The stratification is used by hap.py to report per-GC-bin accuracy.

### 2.6 Advanced Stratification and ClinVar BEDs

```bash
# pwd: variant-calling-benchmark

# Prepare GENCODE CDS, GC, and Coverage stratification BEDs
bash evaluation/prepare_stratification_hg38.sh

# Prepare %GC in CDS Region strata (run this after prepare_stratification_hg38.sh)
bash evaluation/prepare_cds_gc_stratification.sh

# Download and filter ClinVar to chr22 pathogenic variants
bash evaluation/prepare_clinvar.sh
```

`chr22.sdf` is required before running hap.py with `--engine vcfeval`.

### 2.4 Sentieon Assets

Shared-BAM DNAscope and the optional raw-FASTQ DNAscope runs require a valid license. The `07_call_dnascope*.sh` and `07_call_dnascope_fastq*.sh` scripts invoke the local `sentieon-cli` command directly, and `sentieon-cli dnascope` shells out to `sentieon driver`, so both `sentieon-cli` and `sentieon` must be available in `PATH`. For `sentieon-cli dnascope`, `DNASCOPE_WGS_MODEL` must point to the downloaded `.bundle` file, not an unpacked directory.

```bash
# pwd: variant-calling-benchmark

# ---- 1. Download Sentieon license & software ----
mkdir -p .tools/sentieon
# Copy license file vào .tools/sentieon/
cp /path/to/Hanoi_University_Of_Science_And_Technology_eval.lic .tools/sentieon/

wget -O .tools/sentieon/sentieon-genomics-202503.02.tar.gz \
  "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-202503.02.tar.gz"
tar -xzf .tools/sentieon/sentieon-genomics-202503.02.tar.gz -C .tools/sentieon/

# ---- 2. Install sentieon-cli (from local source) ----
pip install .tools/sentieon-cli-src/

# ---- 3. Set environment variables ----
export SENTIEON_LICENSE="$PWD/.tools/sentieon/Hanoi_University_Of_Science_And_Technology_eval.lic"
export SENTIEON_BIN_DIR="$PWD/.tools/sentieon/sentieon-genomics-202503.02/bin"
export PATH="${SENTIEON_BIN_DIR}:$PATH"

# Verify both commands are reachable
command -v sentieon-cli
command -v sentieon

# ---- 4. Download DNAscope model bundles ----
mkdir -p data/reference/models
wget -O data/reference/models/SentieonIlluminaWGS2.2.bundle \
  https://s3.amazonaws.com/sentieon-release/other/SentieonIlluminaWGS2.2.bundle

export DNASCOPE_WGS_MODEL="$PWD/data/reference/models/SentieonIlluminaWGS2.2.bundle"
```

> **Note:** The trial license (up to 128 threads) is valid until **April 15th, 2026**.

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

samtools faidx reference/chr22.fa
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

### 3.2 WGS FASTQs

```bash
# pwd: variant-calling-benchmark/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

for COV in 10 20 30 50; do
  art_illumina \
    -ss HS25 \
    -i "${MUTATED_FASTA}" \
    -p \
    -l 150 \
    -f "${COV}" \
    -m 350 \
    -s 50 \
    -rs 42 \
    -o "${SIM_DIR}/${PREFIX}_${COV}x_" \
    -na

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

  echo "FINAL_BAM=${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/bam_path.sh"
done
```

## 5. Variant Calling

```bash
# pwd: variant-calling-benchmark
for COV in 10 20 30 50; do
  bash pipelines/03_call_hc.sh "${COV}"
  bash pipelines/04_call_dv.sh "${COV}"
  bash pipelines/05_call_strelka.sh "${COV}"
  bash pipelines/06_call_freebayes.sh "${COV}"
  bash pipelines/07_call_dnascope.sh "${COV}"
done
```

### 6.1 Main Track — RTG VCFeval and hap.py

Direct RTG vcfeval evaluation (no Docker needed) and hap.py evaluation (Docker required). Both evaluate all 5 callers × 4 coverages.

```bash
# pwd: variant-calling-benchmark
bash evaluation/eval_vcfeval.sh
bash evaluation/eval_happy.sh    # Optional: run hap.py in parallel
```

Outputs:

- Per-caller results: `results/eval/vcfeval/{COV}x/{caller}/` → `fn.vcf.gz`, `fp.vcf.gz`, `tp.vcf.gz`, `tp-baseline.vcf.gz`, `summary.txt`
- Aggregated stats: `results/eval/vcfeval_all_stats.tsv`


### 6.2 Concordance Analysis

Pairwise comparison between callers using RTG vcfeval:

```bash
# pwd: variant-calling-benchmark

# Step 1: Run pairwise RTG vcfeval for all caller pairs
bash evaluation/concordance/run_concordance.sh

# Step 2: Build concordance matrix from results
python evaluation/concordance/concordance_matrix.py \
    results/eval/concordance \
    results/eval/concordance/concordance_matrix.tsv

# Step 3: Extract variant-level concordance for PCA and UpSet plots
bash evaluation/build_variant_table.sh

# Step 4: Generate heatmap visualizations
python visualization/plot_concordance.py \
    results/eval/concordance/concordance_matrix.tsv \
    results/plots
```

## 7. Visualization

### 7.1 Summary Plots

```bash
# pwd: variant-calling-benchmark
Rscript visualization/benchmark_plots.R results/eval/all_stats.tsv results/plots
python visualization/plot_summary.py results/eval/all_stats.tsv results/plots
```

- **Per-coverage:** F1 bar charts, Precision-Recall scatter, TP/FP/FN counts
- **Cross-coverage:**
  - `grouped_bars_{snp,indel}_*.png` — F1/Recall/Precision grouped by coverage
  - `f1_heatmap_*.png` — Callers × Coverage heatmap
  - `snp_vs_indel_*.png` — Side-by-side SNP vs INDEL F1
  - `variant_breakdown_*.png` — TP/FP/FN per variant type across coverages

### 7.2 Coverage Comparison Plots

Highlights how each caller's performance changes across 10x→50x:

```bash
# pwd: variant-calling-benchmark
python visualization/plot_coverage_comparison.py \
    results/eval/vcfeval_all_stats.tsv results/plots
```

Generates:
- `coverage_lines.png` — F1/Precision/Recall line plots vs coverage
- `coverage_delta_bars.png` — Metric improvement from 10x→50x
- `coverage_radar.png` — Radar chart per coverage
- `coverage_fn_fp_area.png` — Stacked FN+FP area chart
- `coverage_f1_heatmap.png` — F1 callers × coverage heatmap

### 7.3 Advanced Analysis (Pairwise, PCA, Stratified, ClinVar)

These R scripts adapt the analysis pipeline from the reference paper:

```bash
# pwd: variant-calling-benchmark

# 1. Pairwise Wilcoxon Signed-Rank Test Heatmaps
Rscript visualization/01_pairwise_wilcoxon.R

# 2. PCA Scatterplot of Discordant Variants
Rscript visualization/02_pca_scatterplot.R

# 3. Stratified Performance Line Plots (reads hap.py output)
Rscript visualization/03_stratified_plots.R

# 4. ClinVar Pathogenic Variant Detection Analysis
Rscript visualization/04_clinvar_analysis.R

# 5. UpSet Plots
Rscript visualization/upset-plot.R
```

## 8. AlphaGenome Scoring

Score FN/FP variants with AlphaGenome for functional impact analysis:

```bash
# pwd: variant-calling-benchmark
export ALPHAGENOME_API_KEY="your-key"

# Single caller × single coverage
python batch_alphagenome_fn_fp.py --caller gatk --error-type FN --coverage 30

# All coverages at once
python batch_alphagenome_fn_fp.py --caller deepvariant --error-type FP --coverage all

# All 5 callers × FN + FP
for CALLER in gatk deepvariant strelka2 freebayes dnascope; do
  python batch_alphagenome_fn_fp.py --caller $CALLER --error-type FN --coverage all
  python batch_alphagenome_fn_fp.py --caller $CALLER --error-type FP --coverage all
done
```

Output: `results/alphagenome/{caller}_{COV}x_{FN|FP}_variant_scores_tidy.csv`

## 9. Method Notes

- WGS coverages 10x/20x/30x/50x span the regime from low-depth research to production-grade sequencing, enabling depth-impact analysis.

## 10. Future Work

Planned but not implemented in this repo:

- SnpEff / SnpSift annotation
- dbNSFP / ClinVar integration
- ACMG-style FN/FP risk summaries