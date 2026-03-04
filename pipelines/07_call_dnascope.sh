#!/bin/bash
# Sentieon DNAscope — WGS Mode
# Requires: SENTIEON_LICENSE env var or license file
# DNAscope is a machine-learning based variant caller from Sentieon
# Optional: skip gracefully if Sentieon is not available
# Metrics: runtime + Max RSS

set -e

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/bam_path.sh"

CALLER="dnascope"
OUT_DIR="${VARIANT_DIR}/${CALLER}"
TIMEDIR="${LOG_DIR}/time"
mkdir -p "${OUT_DIR}" "${TIMEDIR}"

METRICS="${LOG_DIR}/benchmark_metrics.tsv"
if [[ ! -f "${METRICS}" ]]; then
    echo -e "Caller\tPipeline\tWallClock_sec\tCPU_percent\tMaxRSS_kB" > "${METRICS}"
fi

log_metrics() {
    local caller="$1" pipeline="$2" timefile="$3"
    local wall cpu rss
    wall=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "$timefile")
    cpu=$(awk -F': ' '/Percent of CPU this job got/{print $2}' "$timefile" | tr -d '%')
    rss=$(awk -F': ' '/Maximum resident set size \(kbytes\)/{print $2}' "$timefile")
    local secs
    secs=$(echo "$wall" | awk -F: '{if(NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}')
    echo -e "${caller}\t${pipeline}\t${secs}\t${cpu}\t${rss}" >> "${METRICS}"
    echo "  [METRICS] ${pipeline}: ${wall} wall, ${cpu}% CPU, MaxRSS=${rss} kB"
}

# --- Check Sentieon license ---
if [[ -z "${SENTIEON_LICENSE}" ]]; then
    echo "WARNING: SENTIEON_LICENSE not set. Skipping DNAscope."
    echo "  To use DNAscope, set SENTIEON_LICENSE to your license server or file."
    echo "  Academic users: request free eval at https://www.sentieon.com/free-trial/"
    exit 0
fi

# Check if sentieon is available (either system-installed or Docker)
USE_DOCKER=false
if ! command -v sentieon &>/dev/null; then
    echo "sentieon not found in PATH, will use Docker image: ${SENTIEON_IMAGE}"
    USE_DOCKER=true
fi

echo "=== Running Sentieon DNAscope (WGS) ==="

if [[ "${USE_DOCKER}" == "true" ]]; then
    ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
    ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}" && pwd)
    ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)
    
    BAM_BASENAME=$(basename "${FINAL_BAM}")
    REF_BASENAME=$(basename "${REF_FASTA}")
    DBSNP_BASENAME=$(basename "${DBSNP}")
    
    /usr/bin/time -v -o "${TIMEDIR}/dnascope.time" bash -c "
        set -e
        # Step 1: DNAscope calling
        docker run --rm \
            -e SENTIEON_LICENSE='${SENTIEON_LICENSE}' \
            -v '${ABS_REF_DIR}:/ref:ro' \
            -v '${ABS_PREPROC_DIR}:/input:ro' \
            -v '${ABS_OUT_DIR}:/output' \
            ${SENTIEON_IMAGE} \
            sentieon driver \
            -r '/ref/${REF_BASENAME}' \
            -i '/input/${BAM_BASENAME}' \
            -t ${THREADS} \
            --algo DNAscope \
            -d '/ref/${DBSNP_BASENAME}' \
            '/output/${PREFIX}_DNASCOPE_TMP.vcf'

        # Step 2: Apply DNAscope ML model
        docker run --rm \
            -e SENTIEON_LICENSE='${SENTIEON_LICENSE}' \
            -v '${ABS_REF_DIR}:/ref:ro' \
            -v '${ABS_OUT_DIR}:/output' \
            ${SENTIEON_IMAGE} \
            sentieon driver \
            -r '/ref/${REF_BASENAME}' \
            -t ${THREADS} \
            --algo DNAModelApply \
            --model '/opt/sentieon/share/dnascope_models/SentieonDNAscopeModel1.1.model' \
            -v '/output/${PREFIX}_DNASCOPE_TMP.vcf' \
            '/output/${PREFIX}_DNASCOPE.vcf'
    "
else
    # System-installed sentieon
    /usr/bin/time -v -o "${TIMEDIR}/dnascope.time" bash -c "
        set -e
        sentieon driver \
            -r '${REF_FASTA}' \
            -i '${FINAL_BAM}' \
            -t ${THREADS} \
            --algo DNAscope \
            -d '${DBSNP}' \
            '${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf'

        sentieon driver \
            -r '${REF_FASTA}' \
            -t ${THREADS} \
            --algo DNAModelApply \
            --model /opt/sentieon/share/dnascope_models/SentieonDNAscopeModel1.1.model \
            -v '${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf' \
            '${OUT_DIR}/${PREFIX}_DNASCOPE.vcf'
    "
fi

log_metrics "dnascope" "DNAscope" "${TIMEDIR}/dnascope.time"

# Clean up temp file
rm -f "${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf" "${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf.idx"

echo ""
echo "DNAscope done: ${OUT_DIR}/${PREFIX}_DNASCOPE.vcf"
echo "Metrics: ${METRICS}"
