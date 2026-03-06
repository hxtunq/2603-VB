#!/bin/bash
# FreeBayes + filter — WES Mode
# Same as 06_call_freebayes.sh but with --targets for exome regions
# Metrics: runtime + Max RSS for freebayes calling only

set -e

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/${COV}x_wes/bam_path.sh"

CALLER="freebayes_wes"
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

# --- FreeBayes WES (measured) ---
echo "=== Running FreeBayes (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/freebayes_wes.time" \
    freebayes \
    --standard-filters \
    --targets "${EXOME_BED}" \
    -f "${REF_FASTA}" \
    "${FINAL_BAM}" \
    > "${OUT_DIR}/${PREFIX}.FB.raw.vcf"

log_metrics "freebayes_wes" "FreeBayes_WES" "${TIMEDIR}/freebayes_wes.time"

# --- LeftAlign and split multi-allelics (not measured) ---
echo "Left-aligning variants..."
gatk LeftAlignAndTrimVariants \
    -V "${OUT_DIR}/${PREFIX}.FB.raw.vcf" \
    -R "${REF_FASTA}" \
    --split-multi-allelics \
    -O "${OUT_DIR}/${PREFIX}.FB.split.vcf"

# --- Quality filter (not measured) ---
echo "Filtering variants..."
gatk VariantFiltration \
    -V "${OUT_DIR}/${PREFIX}.FB.split.vcf" \
    -filter "RPR < 1" --filter-name "RPR1" \
    -filter "RPL < 1" --filter-name "RPL1" \
    -filter "SAF < 1" --filter-name "SAF1" \
    -filter "SAR < 1" --filter-name "SAR1" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "QUAL / AO < 10.0" --filter-name "QUALbyAO10" \
    -O "${OUT_DIR}/${PREFIX}_FB_STANDARD.vcf"

echo ""
echo "FB WES done: ${OUT_DIR}/${PREFIX}_FB_STANDARD.vcf"
echo "Metrics: ${METRICS}"
