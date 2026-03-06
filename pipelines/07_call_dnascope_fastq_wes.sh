#!/bin/bash
# Sentieon DNAscope - WES alignment and calling directly from raw FASTQ.

set -euo pipefail

COV="${1:?Usage: $0 <coverage>  (e.g. 50, 100, 200)}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/_sentieon_common.sh"

CALLER="dnascope_fastq_wes"
TIMEFILE="${LOG_DIR}/${COV}x_wes/time/dnascope_fastq_wes.time"
OUTPUT_VCF=""
READGROUP=""
MODEL_BUNDLE=""
CMD=()

sentieon_set_fastq_paths "${COV}" "wes"
sentieon_skip_if_no_license "DNAscope FASTQ WES"
sentieon_require_dnascope_cli_stack
sentieon_require_file "${FASTQ_R1}"
sentieon_require_file "${FASTQ_R2}"
sentieon_require_reference_fasta "${REF_FASTA}"
sentieon_require_bwa_index "${REF_FASTA}"
sentieon_require_file "${DBSNP}"
sentieon_require_vcf_index "${DBSNP}"
sentieon_require_file "${EXOME_BED}"
MODEL_BUNDLE=$(sentieon_resolve_model_bundle "${DNASCOPE_WES_MODEL}" "true")

sentieon_prepare_layout "${COV_SUFFIX}" "${CALLER}"

READGROUP=$(sentieon_build_readgroup "${COV}" "wes")
OUTPUT_VCF="${OUT_DIR}/SS_${CHR_TO_USE}_DNASCOPE_FASTQ_${COV}x_WES.vcf.gz"

CMD=(
    sentieon-cli dnascope
    -r "${REF_FASTA}"
    --r1_fastq "${FASTQ_R1}"
    --r2_fastq "${FASTQ_R2}"
    --readgroups "${READGROUP}"
    -m "${MODEL_BUNDLE}"
    -d "${DBSNP}"
    -b "${EXOME_BED}"
    -t "${THREADS}"
    --assay WES
    "${OUTPUT_VCF}"
)

if [[ "${PCRFREE:-false}" == "true" ]]; then
    CMD+=(--pcr_free)
fi

echo "=== Running Sentieon DNAscope (WES, raw FASTQ) ==="

/usr/bin/time -v -o "${TIMEFILE}" "${CMD[@]}"

sentieon_log_metrics "dnascope_fastq_wes" "DNAscope_FASTQ_WES" "${TIMEFILE}"

echo
echo "DNAscope FASTQ WES done: ${OUTPUT_VCF}"
echo "Metrics: ${METRICS}"
