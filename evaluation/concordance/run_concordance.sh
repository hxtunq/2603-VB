#!/usr/bin/env bash
# Pairwise concordance between callers using RTG vcfeval.
# For each coverage, compares every pair of callers to measure agreement.
#
# Usage: bash evaluation/concordance/run_concordance.sh [COVERAGES...]
#        Default coverages from config if none given.
#
# Requires: rtg (in PATH or via RTG_SH), bcftools, tabix

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../../config/config.sh"

# Callers and their VCF name patterns (must match eval_happy.sh naming)
CALLERS=(gatk deepvariant strelka2 freebayes dnascope dnascope_fastq)

get_vcf() {
    local caller="$1" cov="$2"
    local variant_dir="${VARIANT_DIR}/${cov}x"
    case "${caller}" in
        gatk)            echo "${variant_dir}/gatk/SS_${CHR_TO_USE}_HC_${cov}x_WGS.vcf.gz" ;;
        deepvariant)     echo "${variant_dir}/deepvariant/SS_${CHR_TO_USE}_DV_${cov}x_WGS.vcf.gz" ;;
        strelka2)        echo "${variant_dir}/strelka2/SS_${CHR_TO_USE}_STRELKA_${cov}x_WGS.vcf.gz" ;;
        freebayes)       echo "${variant_dir}/freebayes/SS_${CHR_TO_USE}_FB_${cov}x_WGS.vcf.gz" ;;
        dnascope)        echo "${variant_dir}/dnascope/SS_${CHR_TO_USE}_DNASCOPE_${cov}x_WGS.vcf.gz" ;;
        dnascope_fastq)  echo "${variant_dir}/dnascope_fastq/SS_${CHR_TO_USE}_DNASCOPE_FASTQ_${cov}x_WGS.vcf.gz" ;;
        *)               return 1 ;;
    esac
}

caller_alias() {
    case "$1" in
        gatk)        echo "HC" ;;
        deepvariant) echo "DV" ;;
        strelka2)    echo "ST" ;;
        freebayes)   echo "FB" ;;
        dnascope)        echo "DS" ;;
        dnascope_fastq)  echo "DSF" ;;
        *)               echo "$1" ;;
    esac
}

# Prepare a PASS-filtered + normalized VCF for concordance
prepare_vcf() {
    local vcf="$1"
    local norm="${vcf%.vcf.gz}.pass_norm.vcf.gz"

    if [[ -f "${norm}" ]]; then
        echo "${norm}"
        return
    fi

    if [[ ! -f "${vcf}" ]]; then
        echo ""
        return
    fi

    bcftools view -f 'PASS,.' "${vcf}" \
        | bcftools norm -m -both -f "${REF_FASTA}" \
        | bgzip -c > "${norm}"
    tabix -f -p vcf "${norm}"
    echo "${norm}"
}

# RTG command — use docker or system rtg
run_rtg_vcfeval() {
    local baseline="$1" calls="$2" outdir="$3" bed="$4"

    if command -v rtg &>/dev/null; then
        rtg vcfeval \
            -b "${baseline}" \
            -c "${calls}" \
            -t "${RTG_SDF}" \
            --bed-regions "${bed}" \
            -o "${outdir}" \
            -m combine \
            --squash-ploidy
    else
        echo "ERROR: 'rtg' not found. Install RTG Tools or add to PATH." >&2
        return 1
    fi
}

COVS=("${@:-${COVERAGES_WGS[@]}}")
if [[ ${#COVS[@]} -eq 0 ]]; then
    COVS=("${COVERAGES_WGS[@]}")
fi

CONCORD_ROOT="${EVAL_DIR}/concordance"
mkdir -p "${CONCORD_ROOT}"

echo "=== Pairwise Concordance Analysis ==="
echo "Callers: ${CALLERS[*]}"
echo "Coverages: ${COVS[*]}"
echo ""

for COV in "${COVS[@]}"; do
    echo "--- Coverage: ${COV}x ---"
    COV_DIR="${CONCORD_ROOT}/${COV}x"
    mkdir -p "${COV_DIR}"

    N=${#CALLERS[@]}
    for (( i=0; i<N; i++ )); do
        for (( j=i+1; j<N; j++ )); do
            A="${CALLERS[$i]}"
            B="${CALLERS[$j]}"
            A_ALIAS=$(caller_alias "$A")
            B_ALIAS=$(caller_alias "$B")

            VCF_A=$(get_vcf "$A" "$COV")
            VCF_B=$(get_vcf "$B" "$COV")

            NORM_A=$(prepare_vcf "${VCF_A}")
            NORM_B=$(prepare_vcf "${VCF_B}")

            if [[ -z "${NORM_A}" || -z "${NORM_B}" ]]; then
                echo "  SKIP: ${A_ALIAS} vs ${B_ALIAS} — VCF not found"
                continue
            fi

            PAIR_DIR="${COV_DIR}/${A_ALIAS}_vs_${B_ALIAS}"
            if [[ -d "${PAIR_DIR}" ]]; then
                echo "  EXISTS: ${A_ALIAS} vs ${B_ALIAS} — skipping"
                continue
            fi

            echo "  Running: ${A_ALIAS} vs ${B_ALIAS}"
            run_rtg_vcfeval "${NORM_A}" "${NORM_B}" "${PAIR_DIR}" "${TRUTH_BED}" || {
                echo "  WARN: ${A_ALIAS} vs ${B_ALIAS} failed"
                rm -rf "${PAIR_DIR}"
                continue
            }
        done
    done
done

echo ""
echo "Concordance runs complete. Results in: ${CONCORD_ROOT}/"
echo "Next: python evaluation/concordance/concordance_matrix.py ${CONCORD_ROOT} ${CONCORD_ROOT}/concordance_matrix.tsv"
