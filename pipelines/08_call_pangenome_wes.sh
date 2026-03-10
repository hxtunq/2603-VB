#!/bin/bash
# Sentieon Pangenome DNAscope - WES mode on the shared dedup BAM.
# Uses pangenome graph alignment (vg giraffe) + DNAscope variant calling.

set -euo pipefail

COV="${1:?Usage: $0 <coverage>  (e.g. 50, 100, 200)}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/_sentieon_common.sh"

BAM_PATH_FILE="${PREPROC_DIR}/${COV}x_wes/bam_path.sh"
sentieon_require_file "${BAM_PATH_FILE}"
source "${BAM_PATH_FILE}"

CALLER="pangenome_wes"
TIMEFILE="${LOG_DIR}/${COV}x_wes/time/pangenome_wes.time"
OUTPUT_VCF=""
PANGENOME_MODEL_BUNDLE=""
CMD=()

sentieon_skip_if_no_license "Pangenome DNAscope WES"
sentieon_require_pangenome_stack
sentieon_require_file "${FINAL_BAM}"
sentieon_require_reference_fasta "${REF_FASTA}"
sentieon_require_bwa_index "${REF_FASTA}"
sentieon_require_file "${DBSNP}"
sentieon_require_vcf_index "${DBSNP}"
sentieon_require_file "${EXOME_BED}"
sentieon_require_file "${PANGENOME_GBZ}"
sentieon_require_file "${PANGENOME_HAPL}"
sentieon_require_file "${PANGENOME_POP_VCF}"
sentieon_require_vcf_index "${PANGENOME_POP_VCF}"
PANGENOME_MODEL_BUNDLE=$(sentieon_resolve_model_bundle "${DNASCOPE_PANGENOME_MODEL}")

sentieon_prepare_layout "${COV}x_wes" "${CALLER}"

OUTPUT_VCF="${OUT_DIR}/SS_${CHR_TO_USE}_PANGENOME_${COV}x_WES.vcf.gz"

CMD=(
    sentieon-cli sentieon-pangenome
    -r "${REF_FASTA}"
    -i "${FINAL_BAM}"
    --gbz "${PANGENOME_GBZ}"
    --hapl "${PANGENOME_HAPL}"
    --pop_vcf "${PANGENOME_POP_VCF}"
    -m "${PANGENOME_MODEL_BUNDLE}"
    -d "${DBSNP}"
    -b "${EXOME_BED}"
    -t "${THREADS}"
    --bam_format
    --skip_pangenome_name_checks
    "${OUTPUT_VCF}"
)

if [[ "${PCRFREE:-false}" == "true" ]]; then
    CMD+=(--pcr_free)
fi

echo "=== Running Sentieon Pangenome DNAscope (WES, shared BAM) ==="

/usr/bin/time -v -o "${TIMEFILE}" "${CMD[@]}"

sentieon_log_metrics "pangenome_wes" "Pangenome_DNAscope_WES" "${TIMEFILE}"

echo
echo "Pangenome DNAscope WES done: ${OUTPUT_VCF}"
echo "Metrics: ${METRICS}"
