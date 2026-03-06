#!/bin/bash
set -euo pipefail
#IMPORTANT NOTE: samtools tview is interactive. This script cannot be fully automated.
# Make conda available in non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"
# ============================================================
# chr22 Depth Benchmarking & Variant Analysis 
# ============================================================

# ----------------------------
# Environment & paths
# ----------------------------
conda activate ngs1
cd ~/giab_project/results

REF=~/giab_project/reference/hg38/Homo_sapiens_assembly38.fasta
EXOME_BED=../data/bed/Exome-Agilent_V5_chr22.bed
BASE_BAM=HG002.hg38.bam

# ============================================================
# Step 1: Calculate baseline exome depth (chr22)
# ============================================================
echo "Calculating baseline exome depth..."
samtools depth -b $EXOME_BED $BASE_BAM | \
  awk '{sum+=$3; n++} END {print "Actual exome depth:", sum/n}'

# ============================================================
# Step 2: Downsample BAM to target depths
# Baseline ~220x → probabilities pre-calculated
# ============================================================
# Downsample BAM to target depths
declare -A PVALS=(
  [2x]=0.009
  [10x]=0.047
  [40x]=0.187
  [80x]=0.374
)

for depth in 2x 10x 40x 80x; do
  echo "Downsampling to $depth..."
  gatk DownsampleSam \
    -I "$BASE_BAM" \
    -O "HG002.chr22.${depth}.bam" \
    -P "${PVALS[$depth]}" \
    --RANDOM_SEED 42 \
    --STRATEGY HighAccuracy
done

# Index downsampled BAMs
echo "Indexing downsampled BAMs..."
for bam in HG002.chr22.*x.bam; do
  samtools index $bam
done
# ============================================================
# Step 3: Verify achieved depth
# ============================================================
#Verifying achieved depths
for bam in HG002.chr22.*x.bam; do
  echo -n "$bam -> "
  samtools depth -b $EXOME_BED $bam | \
    awk '{sum+=$3; n++} END {print "Exome DP:", sum/n}'
done

# ============================================================
# Step 4: Variant Calling (GATK HaplotypeCaller)
# ============================================================
#Running variant calling.
for bam in HG002.chr22.*x.bam; do
  name=${bam%.bam}
  gatk --java-options "-Xmx4G" HaplotypeCaller \
    -R $REF \
    -I $bam \
    -L chr22 \
    -O ${name}.vcf.gz
done

# Index VCF files
echo "Indexing VCF files..."
for vcf in HG002.chr22.*x.vcf.gz; do
  tabix -p vcf $vcf
done
# ============================================================
# Step 5: Variant QC (bcftools stats)
# ============================================================
#Running VCF QC
for vcf in HG002.chr22.*x.vcf.gz; do
  name=${vcf%.vcf.gz}
  bcftools stats $vcf > ${name}.stats.txt
  echo "==== $name QC summary ===="
  grep -E "number of SNPs|number of indels|TSTV" ${name}.stats.txt
  echo ""
done

# ============================================================
# Step 6: Validate VCFs
# ============================================================
#Validating VCF files
for vcf in HG002.chr22.*x.vcf.gz; do
  gatk ValidateVariants \
    -R $REF \
    -V $vcf \
    --warn-on-errors
done

# ============================================================
# Step 7: Compare variants across depths
# ============================================================
#Comparing variants across depths

# False positives (2x only)
bcftools isec -C HG002.chr22.2x.vcf.gz HG002.chr22.80x.vcf.gz -p fp_candidates

# False negatives (missed at low depth)
bcftools isec -C HG002.chr22.80x.vcf.gz HG002.chr22.2x.vcf.gz -p fn_candidates

# Shared variants across all depths
bcftools isec -n=4 HG002.chr22.{2,10,40,80}x.vcf.gz -p shared_all

# ============================================================
# Step 8: Extract representative variants
# ============================================================

# Example false positive (high QUAL, low depth)
FP_POS=17834647
bcftools view -H -r chr22:$FP_POS HG002.chr22.2x.vcf.gz
bcftools view -H -r chr22:$FP_POS HG002.chr22.80x.vcf.gz

# Example false negative (real variant)
FN_POS=15470624
bcftools view -H -r chr22:$FN_POS fn_candidates/0000.vcf

# Example shared variant (confidence increases with depth)
SHARED_POS=15628531
bcftools view -H -r chr22:$SHARED_POS shared_all/0000.vcf

# ============================================================
# Step 9: Visualize variants 
# ============================================================
# NOTE: Positions manually selected from Step 8 head output
# - These represent high-quality but false calls at low coverage
# - Demonstrates how low coverage can produce confident-looking errors
# Example false positive (high QUAL, low depth)
FP_POS=17834647
# False positive visualization
samtools tview -p chr22:$FP_POS HG002.chr22.2x.bam  $REF
samtools tview -p chr22:$FP_POS HG002.chr22.80x.bam $REF

# Example false negative (real variant)
FN_POS=15470624

# False negative visualization
samtools tview -p chr22:$FN_POS HG002.chr22.2x.bam  $REF
samtools tview -p chr22:$FN_POS HG002.chr22.80x.bam $REF

# Example shared variant (confidence increases with depth)
SHARED_POS=15628531
# Shared variant visualization
samtools tview -p chr22:$SHARED_POS HG002.chr22.2x.bam  $REF
samtools tview -p chr22:$SHARED_POS HG002.chr22.80x.bam $REF
# ============================================================
# Step 10: Compare indels across depths
# ============================================================
# Comparing indels across depths

# Extract indels from all depths
bcftools view -v indels HG002.chr22.2x.vcf.gz | bcftools view -Oz -o HG002.chr22.2x.indels.vcf.gz
bcftools view -v indels HG002.chr22.10x.vcf.gz | bcftools view -Oz -o HG002.chr22.10x.indels.vcf.gz
bcftools view -v indels HG002.chr22.40x.vcf.gz | bcftools view -Oz -o HG002.chr22.40x.indels.vcf.gz
bcftools view -v indels HG002.chr22.80x.vcf.gz | bcftools view -Oz -o HG002.chr22.80x.indels.vcf.gz

# Index all indel files
bcftools index HG002.chr22.2x.indels.vcf.gz
bcftools index HG002.chr22.10x.indels.vcf.gz
bcftools index HG002.chr22.40x.indels.vcf.gz
bcftools index HG002.chr22.80x.indels.vcf.gz

# Unique to 2x (not in 80x)
bcftools isec -C HG002.chr22.2x.indels.vcf.gz HG002.chr22.80x.indels.vcf.gz -p unique_2x_indels

# Unique to 80x (not in 2x)
bcftools isec -C HG002.chr22.80x.indels.vcf.gz HG002.chr22.2x.indels.vcf.gz -p unique_80x_indels

# Shared indels across all depths (2x, 10x, 40x, 80x)
bcftools isec -n=4 HG002.chr22.2x.indels.vcf.gz HG002.chr22.10x.indels.vcf.gz HG002.chr22.40x.indels.vcf.gz HG002.chr22.80x.indels.vcf.gz -p shared_indels_all
# ============================================================
# Step 11: View representative indels directly
# ============================================================

# View unique to 2x indels (likely false positives)
echo "=== Unique to 2x indels ==="
bcftools view -H unique_2x_indels/0000.vcf | head -3

# View unique to 80x indels (likely false negatives at 2x)
echo "=== Unique to 80x indels ==="
bcftools view -H unique_80x_indels/0000.vcf | head -3

# View shared indels (present in all depths)
echo "=== Shared indels ==="
bcftools view -H shared_indels_all/0000.vcf | head -3

# Special case: Allelic complexity example
COMPLEX_POS=34997184
echo "=== Allelic complexity at position: $COMPLEX_POS ==="
bcftools view -H -r chr22:$COMPLEX_POS HG002.chr22.2x.indels.vcf.gz
bcftools view -H -r chr22:$COMPLEX_POS HG002.chr22.80x.indels.vcf.gz
# ============================================================
# Step 12: Visualize indels (samtools tview)
# ============================================================
# NOTE: Positions manually selected from Step 11 head output
# - These represent high-quality but false calls at low coverage
# - Demonstrates how low coverage can produce confident-looking errors
# Visualize unique to 2x (FALSE POSITIVE - high QUAL=51.6 but wrong)
UNIQUE_2X_POS=20110042
echo "=== Visualizing unique to 2x indel at chr22:$UNIQUE_2X_POS ==="
samtools tview -p chr22:$UNIQUE_2X_POS HG002.chr22.2x.bam $REF
samtools tview -p chr22:$UNIQUE_2X_POS HG002.chr22.80x.bam $REF

# Visualize unique to 80x (FALSE NEGATIVE - real indel missed at low coverage)
UNIQUE_80X_POS=15572440
echo "=== Visualizing unique to 80x indel at chr22:$UNIQUE_80X_POS ==="
samtools tview -p chr22:$UNIQUE_80X_POS HG002.chr22.2x.bam $REF
samtools tview -p chr22:$UNIQUE_80X_POS HG002.chr22.80x.bam $REF

# Visualize allelic complexity (GENOTYPE ERROR - low coverage masks true complexity)
#COMPLEX_POS=34997184
echo "=== Visualizing allelic complexity at chr22:$COMPLEX_POS ==="
samtools tview -p chr22:$COMPLEX_POS HG002.chr22.2x.bam $REF
samtools tview -p chr22:$COMPLEX_POS HG002.chr22.80x.bam $REF

