#!/bin/bash
# Strelka2 Germline via Docker — with Manta candidate indels
# Based on: caller_benchmark-main/pipelines/call_strelka.sh
# Manta provides candidate indel sites to improve Strelka2 INDEL recall
# Metrics: runtime + Max RSS for Manta+Strelka combined

set -e

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/bam_path.sh"

CALLER="strelka2"
OUT_DIR="${VARIANT_DIR}/${CALLER}"
TIMEDIR="${LOG_DIR}/time"
mkdir -p "${OUT_DIR}"/{mantawd,strelkawd} "${TIMEDIR}"

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

ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}" && pwd)
ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)

BAM_BASENAME=$(basename "${FINAL_BAM}")
REF_BASENAME=$(basename "${REF_FASTA}")

# ==========================================================================
# Step 1: Manta — generate candidate indel sites
# ==========================================================================
echo "=== Running Manta for candidate indels ==="
/usr/bin/time -v -o "${TIMEDIR}/manta.time" bash -c "
    set -e
    # Clean previous manta run if exists
    rm -rf '${ABS_OUT_DIR}/mantawd/workspace' '${ABS_OUT_DIR}/mantawd/results'

    # Configure Manta
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${MANTA_IMAGE} \
        configManta.py \
        --bam '/input/${BAM_BASENAME}' \
        --referenceFasta '/ref/${REF_BASENAME}' \
        --runDir '/output/mantawd'

    # Run Manta
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${MANTA_IMAGE} \
        /output/mantawd/runWorkflow.py -m local -j ${THREADS}
"
log_metrics "strelka2" "Manta" "${TIMEDIR}/manta.time"

# Check for Manta candidate indels
CANDIDATE_INDELS="${OUT_DIR}/mantawd/results/variants/candidateSmallIndels.vcf.gz"
if [[ ! -f "${CANDIDATE_INDELS}" ]]; then
    echo "WARNING: Manta candidateSmallIndels not found, running Strelka2 without it"
    INDEL_CANDIDATES_FLAG=""
else
    echo "Manta candidate indels: ${CANDIDATE_INDELS}"
    INDEL_CANDIDATES_FLAG="--indelCandidates '/output/mantawd/results/variants/candidateSmallIndels.vcf.gz'"
fi

# ==========================================================================
# Step 2: Strelka2 — germline variant calling with Manta indel candidates
# ==========================================================================
echo "=== Running Strelka2 (configure + run) ==="
/usr/bin/time -v -o "${TIMEDIR}/strelka2.time" bash -c "
    set -e
    # Clean previous strelka run if exists
    rm -rf '${ABS_OUT_DIR}/strelkawd/workspace' '${ABS_OUT_DIR}/strelkawd/results'

    # Configure Strelka2 with Manta candidate indels
    docker run --rm \
        -v '${ABS_REF_DIR}:/ref:ro' \
        -v '${ABS_PREPROC_DIR}:/input:ro' \
        -v '${ABS_OUT_DIR}:/output' \
        ${STRELKA2_IMAGE} \
        configureStrelkaGermlineWorkflow.py \
        --bam '/input/${BAM_BASENAME}' \
        --referenceFasta '/ref/${REF_BASENAME}' \
        ${INDEL_CANDIDATES_FLAG} \
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
   "${OUT_DIR}/${PREFIX}_STRELKA_STANDART.vcf.gz"

echo ""
echo "Strelka2 done (with Manta candidate indels):"
echo "  Manta candidates: ${CANDIDATE_INDELS}"
echo "  Strelka output:   ${OUT_DIR}/${PREFIX}_STRELKA_STANDART.vcf.gz"
echo "Metrics: ${METRICS}"
