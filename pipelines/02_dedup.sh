#!/bin/bash
# Sort BAM and Mark Duplicates
# Based on: caller_benchmark-main/pipelines/dedup.sh

set -e

source "$(dirname "$0")/../config/config.sh"

INPUT_BAM="${PREPROC_DIR}/${PREFIX}_BWA.bam"
SORTED_BAM="${PREPROC_DIR}/${PREFIX}_BWA.sorted.bam"
DEDUP_BAM="${PREPROC_DIR}/${PREFIX}_BWA.dedup.bam"

echo "Sorting BAM..."
samtools sort -@ ${THREADS} -o "${SORTED_BAM}" "${INPUT_BAM}"
samtools index "${SORTED_BAM}"

echo "Marking duplicates..."
gatk --java-options "${JAVA_OPTS}" MarkDuplicates \
    -I "${SORTED_BAM}" \
    -O "${DEDUP_BAM}" \
    -M "${PREPROC_DIR}/${PREFIX}.MD.metrics"

samtools index "${DEDUP_BAM}"

# Save BAM path for downstream scripts
echo "export FINAL_BAM=\"${DEDUP_BAM}\"" > "${PREPROC_DIR}/bam_path.sh"

echo "Dedup done: ${DEDUP_BAM}"
