#!/bin/bash
set -euo pipefail

# Make conda available in non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"
# ============================================================
# Variant Calling and Benchmarking at 40× Coverage (chr22)
# ============================================================

PROJECT=~/giab_project
RESULTS=$PROJECT/results
REF=$PROJECT/reference/hg38/Homo_sapiens_assembly38.fasta
EXOME_BED=$PROJECT/data/bed/Exome-Agilent_V5_chr22.bed
TRUTH_VCF=$PROJECT/data/truth_vcf/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
TRUTH_BED=$PROJECT/data/bed/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

cd $RESULTS
echo "Starting variant calling at 40× coverage..."

# ------------------------------------------------------------
# 1) Variant Calling with DeepVariant (Docker)
# ------------------------------------------------------------
echo "Running DeepVariant..."

docker run \
  -v "${HOME}/giab_project":"/data" \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref=/data/reference/hg38/Homo_sapiens_assembly38.fasta \
  --reads=/data/results/HG002.chr22.40x.bam \
  --regions=chr22 \
  --output_vcf=/data/results/HG002.chr22.40x.deepvariant.vcf.gz \
  --num_shards=4
tabix -p vcf HG002.chr22.40x.deepvariant.vcf.gz
# ------------------------------------------------------------
# 2) Variant Calling with FreeBayes
# ------------------------------------------------------------
conda activate ngs1 
echo "Running FreeBayes..."

freebayes \
  -f $REF \
  -t $EXOME_BED \
  --min-coverage 5 \
  --min-base-quality 20 \
  --min-mapping-quality 20 \
  HG002.chr22.40x.bam \
  | bgzip -c > HG002.chr22.40x.freebayes.vcf.gz

# Index FreeBayes VCF
tabix -p vcf HG002.chr22.40x.freebayes.vcf.gz

# ------------------------------------------------------------
# 3) Benchmarking with hap.py (Truth-based evaluation)
# ------------------------------------------------------------
echo "Running hap.py benchmarking..."

conda activate hapy_py27

# DeepVariant benchmarking
hap.py \
  $TRUTH_VCF \
  HG002.chr22.40x.deepvariant.vcf.gz \
  -f $TRUTH_BED \
  -T $EXOME_BED \
  -r $REF \
  -o HG002.chr22.40x.deepvariant.happy \
  --threads 4

# FreeBayes benchmarking
hap.py \
  $TRUTH_VCF \
  HG002.chr22.40x.freebayes.vcf.gz \
  -f $TRUTH_BED \
  -T $EXOME_BED \
  -r $REF \
  -o HG002.chr22.40x.freebayes.happy \
  --threads 4

# ------------------------------------------------------------
# 4) Compile Performance Metrics (Recall / Precision / F1)
# ------------------------------------------------------------
echo "Compiling caller comparison metrics..."

echo "Coverage,Caller,Type,TP,FN,FP,Recall,Precision,F1" > caller_comparison.csv

# GATK
awk -F',' -v caller="GATK" '
($1=="SNP" || $1=="INDEL") && $2=="ALL" {
  print "40x,"caller","$1","$4","$5","$7","$11","$12","$14
}' HG002.chr22.40x.happy.summary.csv >> caller_comparison.csv

# DeepVariant
awk -F',' -v caller="DeepVariant" '
($1=="SNP" || $1=="INDEL") && $2=="ALL" {
  print "40x,"caller","$1","$4","$5","$7","$11","$12","$14
}' HG002.chr22.40x.deepvariant.happy.summary.csv >> caller_comparison.csv

# FreeBayes
awk -F',' -v caller="FreeBayes" '
($1=="SNP" || $1=="INDEL") && $2=="ALL" {
  print "40x,"caller","$1","$4","$5","$7","$11","$12","$14
}' HG002.chr22.40x.freebayes.happy.summary.csv >> caller_comparison.csv

# Display summary
column -t -s',' caller_comparison.csv

# ------------------------------------------------------------
# 5) Ts/Tv Ratio Quality Assessment
# ------------------------------------------------------------
echo "Calculating Ts/Tv ratios..."

echo "Caller,Truth_TsTv,Query_TsTv,Difference" > tstv_comparison.csv

# GATK
awk -F',' '$1=="SNP" && $2=="ALL" {
  print "GATK,"$15","$16","$16-$15
}' HG002.chr22.40x.happy.summary.csv >> tstv_comparison.csv

# DeepVariant
awk -F',' '$1=="SNP" && $2=="ALL" {
  print "DeepVariant,"$15","$16","$16-$15
}' HG002.chr22.40x.deepvariant.happy.summary.csv >> tstv_comparison.csv

# FreeBayes
awk -F',' '$1=="SNP" && $2=="ALL" {
  print "FreeBayes,"$15","$16","$16-$15
}' HG002.chr22.40x.freebayes.happy.summary.csv >> tstv_comparison.csv

column -t -s',' tstv_comparison.csv

echo "Variant calling and benchmarking completed successfully."
