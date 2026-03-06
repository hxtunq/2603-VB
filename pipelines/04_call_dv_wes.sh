#!/bin/bash
# DeepVariant via Docker — WES Mode
# Based on: 04_call_dv.sh but with --model_type=WES and --regions
# Metrics: runtime + Max RSS for docker run

set -e

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/${COV}x_wes/bam_path.sh"

CALLER="deepvariant_wes"
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

# Validate exome BED
if [[ ! -f "${EXOME_BED}" ]]; then
    echo "ERROR: EXOME_BED not found: ${EXOME_BED}"
    exit 1
fi

ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}/${COV}x_wes" && pwd)
ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)

BAM_BASENAME=$(basename "${FINAL_BAM}")
REF_BASENAME=$(basename "${REF_FASTA}")
BED_BASENAME=$(basename "${EXOME_BED}")
OUTPUT_PREFIX="SS_${CHR_TO_USE}_DV_${COV}x_WES"

echo "=== Running DeepVariant ${DEEPVARIANT_VERSION} (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/deepvariant_wes.time" \
    docker run --rm \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_PREPROC_DIR}:/input:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    ${DEEPVARIANT_IMAGE} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref="/ref/${REF_BASENAME}" \
    --reads="/input/${BAM_BASENAME}" \
    --regions="/ref/${BED_BASENAME}" \
    --output_vcf="/output/${OUTPUT_PREFIX}.vcf" \
    --output_gvcf="/output/${OUTPUT_PREFIX}.g.vcf" \
    --num_shards=${THREADS}

log_metrics "deepvariant_wes" "DeepVariant_WES" "${TIMEDIR}/deepvariant_wes.time"

echo ""
echo "DV WES done: ${OUT_DIR}/${OUTPUT_PREFIX}.vcf"
echo "Metrics: ${METRICS}"
