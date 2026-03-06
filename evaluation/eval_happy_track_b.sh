#!/usr/bin/env bash
# Evaluate Track B callsets with hap.py — supplementary DNAscope pangenome results.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/_happy_common.sh"

WGS_CALLERS="dnascope_pangenome"
WES_CALLERS="dnascope_pangenome_wes"

get_vcf() {
    local caller="$1"
    local variant_dir="$2"

    case "${caller}" in
        dnascope_pangenome)     printf '%s\n' "${variant_dir}/dnascope_pangenome/${PREFIX}_DNASCOPE_PANGENOME.vcf.gz" ;;
        dnascope_pangenome_wes) printf '%s\n' "${variant_dir}/dnascope_pangenome_wes/${PREFIX}_DNASCOPE_PANGENOME.vcf.gz" ;;
        *)                      return 1 ;;
    esac
}

happy_require_command docker
happy_require_command bgzip
happy_require_command tabix
happy_require_file "${TRUTH_VCF}"
happy_require_file "${REF_FASTA}"
happy_require_file "${TRUTH_BED}"
happy_require_file "${RTG_SDF}"

mkdir -p "${EVAL_TRACK_B_DIR}"

for COV in "${COVERAGES_WGS[@]}"; do
    VARIANT_DIR_COV="${VARIANT_DIR}/${COV}x"

    for caller in ${WGS_CALLERS}; do
        VCF=$(get_vcf "${caller}" "${VARIANT_DIR_COV}")
        if [[ -f "${VCF}" ]] || [[ -f "${VCF}.gz" ]]; then
            happy_make_comparison "${VCF}" "${caller}" "${COV}" "WGS" "${EVAL_TRACK_B_DIR}"
        else
            echo "SKIP: ${caller} @ ${COV}x WGS - VCF not found: ${VCF}"
        fi
    done
done

for COV in "${COVERAGES_WES[@]}"; do
    VARIANT_DIR_COV="${VARIANT_DIR}/${COV}x_wes"

    for caller in ${WES_CALLERS}; do
        VCF=$(get_vcf "${caller}" "${VARIANT_DIR_COV}")
        if [[ -f "${VCF}" ]] || [[ -f "${VCF}.gz" ]]; then
            happy_make_comparison "${VCF}" "${caller}" "${COV}" "WES" "${EVAL_TRACK_B_DIR}"
        else
            echo "SKIP: ${caller} @ ${COV}x WES - VCF not found: ${VCF}"
        fi
    done
done

bash "${SCRIPT_DIR}/gather_stats.sh" "${EVAL_TRACK_B_DIR}" "${EVAL_TRACK_B_DIR}/all_stats.tsv"

echo
echo "Track B evaluation complete."
echo "Results in: ${EVAL_TRACK_B_DIR}/"
