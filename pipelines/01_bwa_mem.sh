#!/bin/bash
# BWA MEM alignment for simulated data
# Based on: caller_benchmark-main/pipelines/run_bwamem.sh

set -e

source "$(dirname "$0")/../config/config.sh"

DATENOW=$( date )
echo "Alignment started at: $DATENOW"

# Input simulated FASTQ
R1="${SIM_DIR}/${PREFIX}_R1.fq.gz"
R2="${SIM_DIR}/${PREFIX}_R2.fq.gz"

mkdir -p "${PREPROC_DIR}" "${LOG_DIR}"

echo "Aligning ${PREFIX} with BWA MEM..."
bwa mem -M -t ${THREADS} \
    -R "${READ_GROUP}" \
    "${REF_FASTA}" "${R1}" "${R2}" \
    2> "${LOG_DIR}/${PREFIX}.bwa.log" \
    | samtools view -bS - > "${PREPROC_DIR}/${PREFIX}_BWA.bam"

DATENOW=$( date )
echo "Alignment finished at: $DATENOW"
