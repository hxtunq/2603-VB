#!/usr/bin/env bash
# Evaluate Track A callsets with hap.py — shared-BAM comparison only.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/_happy_common.sh"

WGS_CALLERS="gatk deepvariant strelka2 freebayes dnascope"
WES_CALLERS="gatk_wes deepvariant_wes strelka2_wes freebayes_wes dnascope_wes"

get_vcf() {
    local caller="$1"
    local variant_dir="$2"
    local cov="$3"
    local mode="$4"

    case "${caller}" in
        gatk)            printf '%s\n' "${variant_dir}/gatk/SS_${CHR_TO_USE}_HC_${cov}x_${mode}.vcf.gz" ;;
        deepvariant)     printf '%s\n' "${variant_dir}/deepvariant/SS_${CHR_TO_USE}_DV_${cov}x_${mode}.vcf.gz" ;;
        strelka2)        printf '%s\n' "${variant_dir}/strelka2/SS_${CHR_TO_USE}_STRELKA_${cov}x_${mode}.vcf.gz" ;;
        freebayes)       printf '%s\n' "${variant_dir}/freebayes/SS_${CHR_TO_USE}_FB_${cov}x_${mode}.vcf.gz" ;;
        dnascope)        printf '%s\n' "${variant_dir}/dnascope/SS_${CHR_TO_USE}_DNASCOPE_${cov}x_${mode}.vcf.gz" ;;
        gatk_wes)        printf '%s\n' "${variant_dir}/gatk_wes/SS_${CHR_TO_USE}_HC_${cov}x_${mode}.vcf.gz" ;;
        deepvariant_wes) printf '%s\n' "${variant_dir}/deepvariant_wes/SS_${CHR_TO_USE}_DV_${cov}x_${mode}.vcf.gz" ;;
        strelka2_wes)    printf '%s\n' "${variant_dir}/strelka2_wes/SS_${CHR_TO_USE}_STRELKA_${cov}x_${mode}.vcf.gz" ;;
        freebayes_wes)   printf '%s\n' "${variant_dir}/freebayes_wes/SS_${CHR_TO_USE}_FB_${cov}x_${mode}.vcf.gz" ;;
        dnascope_wes)    printf '%s\n' "${variant_dir}/dnascope_wes/SS_${CHR_TO_USE}_DNASCOPE_${cov}x_${mode}.vcf.gz" ;;
        *)               return 1 ;;
    esac
}

happy_require_command docker
happy_require_command bgzip
happy_require_command tabix
happy_require_file "${TRUTH_VCF}"
happy_require_file "${REF_FASTA}"
happy_require_file "${TRUTH_BED}"
happy_require_file "${RTG_SDF}"

mkdir -p "${EVAL_DIR}"

for COV in "${COVERAGES_WGS[@]}"; do
    VARIANT_DIR_COV="${VARIANT_DIR}/${COV}x"

    for caller in ${WGS_CALLERS}; do
        VCF=$(get_vcf "${caller}" "${VARIANT_DIR_COV}" "${COV}" "WGS")
        if [[ -f "${VCF}" ]] || [[ -f "${VCF}.gz" ]]; then
            happy_make_comparison "${VCF}" "${caller}" "${COV}" "WGS" "${EVAL_DIR}"
        else
            echo "SKIP: ${caller} @ ${COV}x WGS - VCF not found: ${VCF}"
        fi
    done
done

for COV in "${COVERAGES_WES[@]}"; do
    VARIANT_DIR_COV="${VARIANT_DIR}/${COV}x_wes"

    for caller in ${WES_CALLERS}; do
        VCF=$(get_vcf "${caller}" "${VARIANT_DIR_COV}" "${COV}" "WES")
        if [[ -f "${VCF}" ]] || [[ -f "${VCF}.gz" ]]; then
            happy_make_comparison "${VCF}" "${caller}" "${COV}" "WES" "${EVAL_DIR}"
        else
            echo "SKIP: ${caller} @ ${COV}x WES - VCF not found: ${VCF}"
        fi
    done
done

bash "${SCRIPT_DIR}/gather_stats.sh" "${EVAL_DIR}" "${EVAL_DIR}/all_stats.tsv"

echo
echo "Track A evaluation complete."
echo "Results in: ${EVAL_DIR}/"
