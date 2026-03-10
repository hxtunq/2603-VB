#!/bin/bash
# Strelka2 Germline via Docker — WES mode (standalone, no Manta)
# Same as 05_call_strelka.sh but with --exome flag and --callRegions

set -e

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/${COV}x_wes/bam_path.sh"

CALLER="strelka2_wes"
OUT_DIR="${VARIANT_DIR}/${COV}x_wes/${CALLER}"
TIMEDIR="${LOG_DIR}/${COV}x_wes/time"
mkdir -p "${OUT_DIR}/strelkawd" "${TIMEDIR}"

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
FINAL_VCF="${OUT_DIR}/SS_${CHR_TO_USE}_STRELKA_${COV}x_WES.vcf.gz"

# Strelka2 requires bgzipped+tabixed BED for --callRegions
CALL_REGIONS_BED="${OUT_DIR}/exome_targets.bed.gz"
if [[ ! -f "${CALL_REGIONS_BED}" ]]; then
    bgzip -c "${EXOME_BED}" > "${CALL_REGIONS_BED}"
    tabix -p bed "${CALL_REGIONS_BED}"
fi
CALL_REGIONS_BASENAME=$(basename "${CALL_REGIONS_BED}")

# ==========================================================================
# Strelka2 — WES mode with --exome (standalone)
# ==========================================================================
echo "=== Running Strelka2 (WES + exome mode) ==="
/usr/bin/time -v -o "${TIMEDIR}/strelka2_wes.time" bash -c "
    set -e
    rm -rf '${ABS_OUT_DIR}/strelkawd/workspace' '${ABS_OUT_DIR}/strelkawd/results'

    # Configure Strelka2 with --exome and --callRegions
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${STRELKA2_IMAGE} \
        configureStrelkaGermlineWorkflow.py \
        --bam '/input/${BAM_BASENAME}' \
        --referenceFasta '/ref/${REF_BASENAME}' \
        --exome \
        --callRegions '/output/${CALL_REGIONS_BASENAME}' \
        --runDir '/output/strelkawd'

    # Run Strelka2
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${STRELKA2_IMAGE} \
        /output/strelkawd/runWorkflow.py -m local -j ${THREADS}
"
log_metrics "strelka2_wes" "Strelka2_WES" "${TIMEDIR}/strelka2_wes.time"

# --- Copy output ---
cp "${OUT_DIR}/strelkawd/results/variants/variants.vcf.gz" \
   "${FINAL_VCF}"

echo ""
echo "Strelka2 WES done (standalone, no Manta):"
echo "  Output: ${FINAL_VCF}"
echo "Metrics: ${METRICS}"
