#!/bin/bash
# Strelka2 Germline via Docker — WGS mode (standalone, no Manta)
# Based on: caller_benchmark-main/pipelines/call_strelka.sh

set -e

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/${COV}x/bam_path.sh"

CALLER="strelka2"
OUT_DIR="${VARIANT_DIR}/${COV}x/${CALLER}"
TIMEDIR="${LOG_DIR}/${COV}x/time"
mkdir -p "${OUT_DIR}/strelkawd" "${TIMEDIR}"

METRICS="${LOG_DIR}/${COV}x/benchmark_metrics.tsv"
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

ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}/${COV}x" && pwd)
ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)

BAM_BASENAME=$(basename "${FINAL_BAM}")
REF_BASENAME=$(basename "${REF_FASTA}")
FINAL_VCF="${OUT_DIR}/SS_${CHR_TO_USE}_STRELKA_${COV}x_WGS.vcf.gz"

# ==========================================================================
# Strelka2 — germline variant calling (standalone)
# ==========================================================================
echo "=== Running Strelka2 (configure + run) ==="
/usr/bin/time -v -o "${TIMEDIR}/strelka2.time" bash -c "
    set -e
    # Clean previous strelka run if exists
    rm -rf '${ABS_OUT_DIR}/strelkawd/workspace' '${ABS_OUT_DIR}/strelkawd/results'

    # Configure Strelka2
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${STRELKA2_IMAGE} \
        configureStrelkaGermlineWorkflow.py \
        --bam '/input/${BAM_BASENAME}' \
        --referenceFasta '/ref/${REF_BASENAME}' \
        --runDir '/output/strelkawd'

    # Run Strelka2
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${STRELKA2_IMAGE} \
        /output/strelkawd/runWorkflow.py -m local -j ${THREADS}
"
log_metrics "strelka2" "Strelka2" "${TIMEDIR}/strelka2.time"

# --- Copy output ---
cp "${OUT_DIR}/strelkawd/results/variants/variants.vcf.gz" \
   "${FINAL_VCF}"

echo ""
echo "Strelka2 WGS done (standalone, no Manta):"
echo "  Output: ${FINAL_VCF}"
echo "Metrics: ${METRICS}"
