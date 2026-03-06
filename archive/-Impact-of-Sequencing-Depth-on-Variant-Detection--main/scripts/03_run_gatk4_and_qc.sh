#!/bin/bash
set -euo pipefail

# Make conda available in non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"
# ------------------------------------------------------------
# Run GATK Cromwell Workflow (REQUIRES cromwell ENV)
# ------------------------------------------------------------
cd ~
conda activate cromwell
docker version

java -Xmx2g -Dconfig.file=cromwell.conf -jar cromwell.jar run \
  ~/giab_project/scripts/processing-for-variant-discovery-gatk4.wdl \
  --inputs ~/giab_project/scripts/processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json

# ------------------------------------------------------------
# Copy Cromwell outputs 
# ------------------------------------------------------------
# NOTE: Cromwell generates a unique workflow ID for each run.
#   Multiple runs create separate output directories with different IDs.
#   This command finds the first matching output file, which may not be
#   the most recent if multiple workflow runs exist. Consider sorting by
#   modification time (-printf '%T@ %p\n' | sort -rn) for production use.
WORKFLOW_DIR=$(find ~/cromwell-executions/PreProcessingForVariantDiscovery_GATK4/ \
  -name "HG002.hg38.bam" -type f | head -1 | xargs dirname)

if [ -z "$WORKFLOW_DIR" ] || [ ! -d "$WORKFLOW_DIR" ]; then
    echo "ERROR: Workflow output directory not found"
    exit 1
fi

echo "Found outputs in: $WORKFLOW_DIR"

# Verify main files exist
for file in HG002.hg38.bam HG002.hg38.bai HG002.hg38.bam.md5; do
    if [ ! -f "$WORKFLOW_DIR/$file" ]; then
        echo "ERROR: $file not found in workflow output"
        exit 1
    fi
done

# Copy main BAM files
cp $WORKFLOW_DIR/HG002.hg38.bam ~/giab_project/results/HG002.hg38.bam
cp $WORKFLOW_DIR/HG002.hg38.bai ~/giab_project/results/HG002.hg38.bai
cp $WORKFLOW_DIR/HG002.hg38.bam.md5 ~/giab_project/results/HG002.hg38.bam.md5

# Copy metrics and recalibration files
cp "$WORKFLOW_DIR/../call-MarkDuplicates/execution/HG002.hg38.duplicate_metrics" \
   ~/giab_project/results/HG002.hg38.duplicate_metrics
cp "$WORKFLOW_DIR/../call-GatherBqsrReports/execution/HG002.hg38.recal_data.csv" \
   ~/giab_project/results/HG002.hg38.recal_data.csv

ls -lh ~/giab_project/results/
# ------------------------------------------------------------
# Basic BAM Statistics 
# ------------------------------------------------------------
cd ~/giab_project/results
conda activate ngs1

ls -lh HG002.hg38.bam*

samtools flagstat HG002.hg38.bam

samtools depth -r chr22 HG002.hg38.bam \
| awk '{sum+=$3; n++} END {print sum/n}'

# ------------------------------------------------------------
# Verify Chr22 Only 
# ------------------------------------------------------------
samtools idxstats HG002.hg38.bam

# ------------------------------------------------------------
# Check Duplicate Rate
# ------------------------------------------------------------
cat HG002.hg38.duplicate_metrics

grep "PERCENT_DUPLICATION" -A 1 HG002.hg38.duplicate_metrics

# ------------------------------------------------------------
# Check BQSR Report 
# ------------------------------------------------------------
head -50 HG002.hg38.recal_data.csv

# ------------------------------------------------------------
# Validate BAM Integrity 
# ------------------------------------------------------------
gatk ValidateSamFile \
  -I HG002.hg38.bam \
  -MODE SUMMARY

# ------------------------------------------------------------
# Coverage Distribution 
# ------------------------------------------------------------
samtools depth -r chr22 HG002.hg38.bam | \
  awk '{sum+=$3; count++} END {
    print "Mean coverage: " sum/count "x"
    print "Total bases covered: " count
  }'

samtools depth -r chr22 HG002.hg38.bam | \
  awk '{print $3}' | \
  sort -n | \
  uniq -c | \
  awk '{print $2"\t"$1}' | \
  head -20

# ------------------------------------------------------------
# Visual Inspection 
# ------------------------------------------------------------
samtools view HG002.hg38.bam chr22:16050000-16051000 | head -5

samtools view -H HG002.hg38.bam | grep "^@RG"

echo "Pipeline finished successfully."

---


