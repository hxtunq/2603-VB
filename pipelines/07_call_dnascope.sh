#!/bin/bash
# Sentieon DNAscope - WGS mode on the shared dedup BAM via sentieon-cli.

set -euo pipefail

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50)}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/_sentieon_common.sh"

BAM_PATH_FILE="${PREPROC_DIR}/${COV}x/bam_path.sh"
sentieon_require_file "${BAM_PATH_FILE}"
source "${BAM_PATH_FILE}"

CALLER="dnascope"
TIMEFILE="${LOG_DIR}/${COV}x/time/dnascope.time"
OUTPUT_VCF=""
MODEL_BUNDLE=""
CMD=()

sentieon_skip_if_no_license "DNAscope"
sentieon_require_dnascope_cli_stack
sentieon_require_file "${FINAL_BAM}"
sentieon_require_reference_fasta "${REF_FASTA}"
sentieon_require_file "${DBSNP}"
sentieon_require_vcf_index "${DBSNP}"
MODEL_BUNDLE=$(sentieon_resolve_model_bundle "${DNASCOPE_WGS_MODEL}")

sentieon_prepare_layout "${COV}x" "${CALLER}"

OUTPUT_VCF="${OUT_DIR}/SS_${CHR_TO_USE}_DNASCOPE_${COV}x_WGS.vcf.gz"

CMD=(
    sentieon-cli dnascope
    -r "${REF_FASTA}"
    -i "${FINAL_BAM}"
    -m "${MODEL_BUNDLE}"
    -d "${DBSNP}"
    -t "${THREADS}"
    --duplicate_marking none
    --assay WGS
    "${OUTPUT_VCF}"
)

if [[ "${PCRFREE:-false}" == "true" ]]; then
    CMD+=(--pcr_free)
fi

echo "=== Running Sentieon DNAscope (WGS, shared BAM) ==="

/usr/bin/time -v -o "${TIMEFILE}" "${CMD[@]}"

sentieon_log_metrics "dnascope" "DNAscope" "${TIMEFILE}"

echo
echo "DNAscope done: ${OUTPUT_VCF}"
echo "Metrics: ${METRICS}"
