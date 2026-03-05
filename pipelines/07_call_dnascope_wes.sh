#!/bin/bash
# Sentieon DNAscope — WES Mode
# Same as 07_call_dnascope.sh but with --interval for exome targets
# Requires: SENTIEON_LICENSE env var

set -e

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/${COV}x_wes/bam_path.sh"

CALLER="dnascope_wes"
OUT_DIR="${VARIANT_DIR}/${COV}x_wes/${CALLER}"
TIMEDIR="${LOG_DIR}/${COV}x_wes/time"
mkdir -p "${OUT_DIR}" "${TIMEDIR}"

METRICS="${LOG_DIR}/${COV}x_wes/benchmark_metrics.tsv"
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

# --- Check ---
if [[ -z "${SENTIEON_LICENSE}" ]]; then
    echo "WARNING: SENTIEON_LICENSE not set. Skipping DNAscope WES."
    exit 0
fi

if [[ ! -f "${EXOME_BED}" ]]; then
    echo "ERROR: EXOME_BED not found: ${EXOME_BED}"
    exit 1
fi

USE_DOCKER=false
if ! command -v sentieon &>/dev/null; then
    USE_DOCKER=true
fi

echo "=== Running Sentieon DNAscope (WES) ==="

if [[ "${USE_DOCKER}" == "true" ]]; then
    ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
    ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}/${COV}x_wes" && pwd)
    ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)
    
    BAM_BASENAME=$(basename "${FINAL_BAM}")
    REF_BASENAME=$(basename "${REF_FASTA}")
    DBSNP_BASENAME=$(basename "${DBSNP}")
    BED_BASENAME=$(basename "${EXOME_BED}")
    
    /usr/bin/time -v -o "${TIMEDIR}/dnascope_wes.time" bash -c "
        set -e
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
            --interval '/ref/${BED_BASENAME}' \
            --algo DNAscope \
            -d '/ref/${DBSNP_BASENAME}' \
            '/output/${PREFIX}_DNASCOPE_TMP.vcf'

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
    /usr/bin/time -v -o "${TIMEDIR}/dnascope_wes.time" bash -c "
        set -e
        sentieon driver \
            -r '${REF_FASTA}' \
            -i '${FINAL_BAM}' \
            -t ${THREADS} \
            --interval '${EXOME_BED}' \
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

log_metrics "dnascope_wes" "DNAscope_WES" "${TIMEDIR}/dnascope_wes.time"

rm -f "${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf" "${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf.idx"

echo ""
echo "DNAscope WES done: ${OUT_DIR}/${PREFIX}_DNASCOPE.vcf"
echo "Metrics: ${METRICS}"
