#!/bin/bash
set -euo pipefail
# This script is designed and validated for Windows Subsystem for Linux (WSL).
# Make conda available in non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"
# ==============================
# GIAB HG002 Chr22 Exome Analysis
# Complete Download Script
# ==============================

PROJECT=~/giab_project
mkdir -p $PROJECT
cd $PROJECT

# ==============================
# 1) Create project structure
# ==============================
echo "[1] Creating project structure"
mkdir -p data/bam data/bed data/truth_vcf data/known_sites
mkdir -p reference/hg38
mkdir -p scripts results 

# ==============================
# 2) Download reference (hg38)
# ==============================
echo "[2] Downloading hg38 reference genome"
cd $PROJECT/reference/hg38

# Core reference files
gsutil -m cp \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict \
  .

# BWA index files
gsutil -m cp \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa \
  .

# ==============================
# 3) Download GIAB truth set
# ==============================
echo "[3] Downloading GIAB truth VCF and BED files"
GIAB_FTP="ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"

cd $PROJECT/data/truth_vcf
wget -c $GIAB_FTP/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c $GIAB_FTP/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

cd $PROJECT/data/bed
wget -c $GIAB_FTP/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# ==============================
# 4) Download Agilent Exome V5 BED
# ==============================
echo "[4] Downloading Agilent V5 exome capture regions"
cd $PROJECT/data/bed
wget -c https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/Exome-Agilent_V5.bed
# Subset to chr22 only
grep -P "^chr22\t" Exome-Agilent_V5.bed > Exome-Agilent_V5_chr22.bed
# ==============================
# 5) Download exome BAM
# ==============================
 echo "[5] Downloading HG002 exome BAM file"
GIAB_BAM_FTP="ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome"

cd $PROJECT/data/bam
wget -c $GIAB_BAM_FTP/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
wget -c $GIAB_BAM_FTP/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
# Rename to simpler names
mv 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam HG002_Oslo_exome.bam
mv 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai HG002_Oslo_exome.bam.bai

# ==============================
# 6) Download known sites
# ==============================
 echo "[6] Downloading known variant sites"
cd $PROJECT/reference/hg38

gsutil -m cp \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
  gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
  .

# ==============================
# 7) Download WDL workflows
# ==============================
 echo "[7] Downloading GATK WDL workflows"
cd $PROJECT/scripts

wget -c https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/6b69dde1877de997effedd3388aeaca929db04da/bam-to-unmapped-bams.wdl
wget -c https://raw.githubusercontent.com/gatk-workflows/seq-format-conversion/6b69dde1877de997effedd3388aeaca929db04da/bam-to-unmapped-bams.inputs.json
wget -c https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/c44603c464fe3cb7d9b82da2a95f844fdeb20e3c/processing-for-variant-discovery-gatk4.wdl
wget -c https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/c44603c464fe3cb7d9b82da2a95f844fdeb20e3c/processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json

# ==============================
# Verify downloads
# ==============================

echo " Verifying critical files "
cd $PROJECT

# Check reference
ls -lh reference/hg38/Homo_sapiens_assembly38.fasta reference/hg38/Homo_sapiens_assembly38.fasta.64.* | head -8

# Check GIAB truth
ls -lh data/truth_vcf/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz*

# Check BED files
ls -lh data/bed/*.bed

# Check BAM
ls -lh data/bam/*.bam*

# Check known sites
ls -lh reference/hg38/*.vcf.gz | head -6

# Check WDL
ls -lh scripts/*.wdl scripts/*.json

echo "Verification Finished"