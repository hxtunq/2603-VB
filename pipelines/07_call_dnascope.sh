#!/bin/bash
# Sentieon DNAscope — WGS mode on the shared dedup BAM.

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

sentieon_skip_if_no_license "DNAscope"
sentieon_require_command docker
sentieon_require_file "${FINAL_BAM}"
sentieon_require_file "${REF_FASTA}"
sentieon_require_file "${DBSNP}"
sentieon_require_file "${DNASCOPE_WGS_MODEL}/dnascope.model"

sentieon_prepare_layout "${COV}x" "${CALLER}"
sentieon_init_license_args
sentieon_set_pcr_arg

ABS_REF_DIR=$(sentieon_abs_dir "${REF_DIR}")
ABS_PREPROC_DIR=$(sentieon_abs_dir "$(dirname "${FINAL_BAM}")")
ABS_OUT_DIR=$(sentieon_abs_dir "${OUT_DIR}")
ABS_MODEL_DIR=$(sentieon_abs_dir "${DNASCOPE_WGS_MODEL}")

REF_BASENAME=$(basename "${REF_FASTA}")
DBSNP_BASENAME=$(basename "${DBSNP}")
BAM_BASENAME=$(basename "${FINAL_BAM}")

echo "=== Running Sentieon DNAscope (WGS, shared BAM) ==="

/usr/bin/time -v -o "${TIMEFILE}" \
    docker run --rm \
    "${SENTIEON_DOCKER_LICENSE_ARGS[@]}" \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_PREPROC_DIR}:/input:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    -v "${ABS_MODEL_DIR}:/model:ro" \
    "${SENTIEON_IMAGE}" \
    bash -lc "
        set -euo pipefail
        sentieon driver \
            -r /ref/${REF_BASENAME} \
            -t ${THREADS} \
            -i /input/${BAM_BASENAME} \
            --algo DNAscope ${SENTIEON_PCR_ARG} \
            --model /model/dnascope.model \
            -d /ref/${DBSNP_BASENAME} \
            /output/${PREFIX}_DNASCOPE_TMP.vcf

        sentieon driver \
            -r /ref/${REF_BASENAME} \
            -t ${THREADS} \
            --algo DNAModelApply \
            --model /model/dnascope.model \
            -v /output/${PREFIX}_DNASCOPE_TMP.vcf \
            /output/${PREFIX}_DNASCOPE.vcf
    "

sentieon_log_metrics "dnascope" "DNAscope" "${TIMEFILE}"

rm -f "${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf" "${OUT_DIR}/${PREFIX}_DNASCOPE_TMP.vcf.idx"

echo
echo "DNAscope done: ${OUT_DIR}/${PREFIX}_DNASCOPE.vcf"
echo "Metrics: ${METRICS}"
