#!/bin/bash
#===============================================================================
# Preprocessing pipeline: FASTQ → analysis-ready BAM
# BWA-MEM → sort → MarkDuplicates → coverage stats
#
# NOTE: BQSR is NOT done here — it's GATK HC-specific (per Barbitoff et al.)
#       Other callers (DeepVariant, Strelka2, FreeBayes) use the dedup BAM directly.
#
# Usage:
#   bash pipelines/00_preprocess.sh <coverage>
#   Example: bash pipelines/00_preprocess.sh 30
#
# Input:  data/simulated/SIMULATED_SAMPLE_chr22_{COV}x_R1.fastq.gz
# Output: results/preprocessing/{COV}x/SIMULATED_SAMPLE_chr22_dedup.bam
#===============================================================================

set -e

source "$(dirname "$0")/../config/config.sh"

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

echo "============================================"
echo "Preprocessing: ${COV}x coverage"
echo "============================================"

# --- Paths ---
R1="${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq.gz"
R2="${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq.gz"
OUTDIR="${RESULTS_DIR}/preprocessing/${COV}x"
LOGDIR="${LOG_DIR}/${COV}x"
mkdir -p "${OUTDIR}" "${LOGDIR}"

if [[ ! -f "${R1}" ]]; then
    echo "ERROR: FASTQ not found: ${R1}"
    exit 1
fi

# ==========================================================================
# 1. FastQC
# ==========================================================================
echo "[1/4] FastQC..."
mkdir -p "${OUTDIR}/fastqc"
fastqc -t ${THREADS} -o "${OUTDIR}/fastqc" "${R1}" "${R2}" \
    2>&1 | tee "${LOGDIR}/fastqc.log"

# ==========================================================================
# 2. BWA-MEM → sorted BAM
# ==========================================================================
echo "[2/4] BWA-MEM alignment + sort..."
bwa mem -t ${THREADS} -M \
    -R "@RG\tID:${PREFIX}_${COV}x\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF_FASTA}" "${R1}" "${R2}" \
    2> "${LOGDIR}/bwa_mem.log" | \
    samtools sort -@ ${THREADS} -m 2G \
    -o "${OUTDIR}/${PREFIX}_aligned.bam" -

samtools index "${OUTDIR}/${PREFIX}_aligned.bam"

# ==========================================================================
# 3. MarkDuplicates
# ==========================================================================
echo "[3/4] MarkDuplicates..."
gatk --java-options "${JAVA_OPTS}" MarkDuplicates \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_dedup.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true \
    2>&1 | tee "${LOGDIR}/markduplicates.log"

# ==========================================================================
# 4. Coverage stats (samtools + mosdepth)
# ==========================================================================
echo "[4/4] Coverage stats..."
samtools stats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
samtools flagstat "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
samtools idxstats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_idxstats.txt"

mosdepth -t ${THREADS} --by 1000 \
    "${OUTDIR}/${PREFIX}_coverage" \
    "${OUTDIR}/${PREFIX}_dedup.bam"

# ==========================================================================
# Export BAM path for downstream caller scripts
# ==========================================================================
echo "FINAL_BAM=${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/bam_path.sh"

echo ""
echo "Done: ${OUTDIR}/${PREFIX}_dedup.bam"
echo "  → This BAM is used by DeepVariant, Strelka2, FreeBayes, DNAscope"
echo "  → GATK HC will run BQSR internally before calling"
