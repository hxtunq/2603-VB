#!/bin/bash
set -euo pipefail
# This script is designed and validated for Windows Subsystem for Linux (WSL).
# Make conda available in non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"
# ============================================================
# BAM → uBAM conversion and chr22 subsetting
# ============================================================

PROJECT=~/giab_project
BAM_DIR=$PROJECT/data/bam

echo "Starting BAM → uBAM conversion..."

# ------------------------------------------------------------
# 1) Run Cromwell (CROMWELL ENV REQUIRED)
# ------------------------------------------------------------
conda activate cromwell
java -jar ~/cromwell.jar run \
  ~/giab_project/scripts/bam-to-unmapped-bams.wdl \
  --inputs ~/giab_project/scripts/bam-to-unmapped-bams.inputs.json
# ------------------------------------------------------------
# 2) Copy uBAM output (explicit path, reproducible)
# ------------------------------------------------------------
# NOTE: Cromwell generates a unique workflow ID for each run.
#   Multiple runs create separate output directories with different IDs.
#   This command finds the first matching output file, which may not be
#   the most recent if multiple workflow runs exist. Consider sorting by
#   modification time (-printf '%T@ %p\n' | sort -rn) for production use.

UBAM=$(find ~/cromwell-executions/BamToUnmappedBams \
  -name "*.unmapped.bam" -type f | head -n 1)

if [ -z "$UBAM" ] || [ ! -f "$UBAM" ]; then
    echo "ERROR: uBAM output not found"
    exit 1
fi

echo "Found uBAM: $UBAM"
cp "$UBAM" $BAM_DIR/HG002.unmapped.bam

# ------------------------------------------------------------
# 3) Validate uBAM (NGS ENV)
# ------------------------------------------------------------
conda activate ngs1

# Check BAM integrity
samtools quickcheck -v HG002.unmapped.bam

# Inspect header (read groups preserved)
samtools view -H HG002.unmapped.bam

# Basic statistics
samtools flagstat HG002.unmapped.bam

# Inspect first reads
samtools view HG002.unmapped.bam | head

# ------------------------------------------------------------
# 4) Extract chr22 read names from aligned BAM
# ------------------------------------------------------------
# Inspect chromosome naming
samtools view -H HG002_Oslo_exome.bam | grep "^@SQ"

# Extract read names mapped to chromosome 22 (hg19 numeric naming)
samtools view HG002_Oslo_exome.bam 22 | \
  cut -f1 | sort -u > chr22_read_names.txt

# Validate extraction
wc -l chr22_read_names.txt
head -5 chr22_read_names.txt
ls -lh chr22_read_names.txt
# ------------------------------------------------------------
# 5) Filter uBAM to chr22 only
# ------------------------------------------------------------
java -Xmx4g -jar ~/picard.jar FilterSamReads \
  I=HG002.unmapped.bam \
  O=HG002_chr22.unmapped.bam \
  READ_LIST_FILE=chr22_read_names.txt \
  FILTER=includeReadList

# ------------------------------------------------------------
# 6) Validate chr22 uBAM
# NOTE:
#   chr22 uBAM is an UNMAPPED BAM by design.
#   - No @SQ lines in header is expected
#   - 0 mapped reads in flagstat is expected
#   Reads were selected using hg19 alignment, but alignment
#   coordinates were intentionally removed to avoid bias.
# ------------------------------------------------------------

samtools quickcheck -v HG002_chr22.unmapped.bam
samtools flagstat HG002_chr22.unmapped.bam
samtools view HG002_chr22.unmapped.bam | head

echo "chr22 uBAM generation completed successfully."
