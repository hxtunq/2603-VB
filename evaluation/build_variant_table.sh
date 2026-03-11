#!/usr/bin/env bash
# =============================================================================
# Build variant-level concordance table from RTG vcfeval outputs.
#
# For each coverage, reads tp.vcf.gz / fp.vcf.gz / fn.vcf.gz per caller and
# produces:
#   1) variant_concordance_{cov}x.tsv   — binary matrix (variant × caller)
#   2) upset_long_callset_{cov}x.csv    — long format FP for UpSet plot
#   3) upset_long_baseline_{cov}x.csv   — long format FN for UpSet plot
#
# Usage:
#   bash evaluation/build_variant_table.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

WGS_CALLERS="${WGS_CALLERS:-gatk deepvariant strelka2 freebayes dnascope}"

caller_alias() {
    case "$1" in
        gatk)        echo "HC" ;;
        deepvariant) echo "DV" ;;
        strelka2)    echo "ST" ;;
        freebayes)   echo "FB" ;;
        dnascope)    echo "DS" ;;
        *)           echo "$1" ;;
    esac
}

# Extract variant IDs from a VCF (CHROM:POS:REF:ALT)
extract_variant_ids() {
    local vcf="$1"
    if [[ ! -f "${vcf}" ]]; then
        return
    fi
    zcat "${vcf}" | awk '!/^#/ {print $1":"$2":"$4":"$5}'
}

VCFEVAL_ROOT="${EVAL_DIR}/vcfeval"
OUT_DIR="${VCFEVAL_ROOT}"

echo "=========================================="
echo " Building Variant Concordance Tables"
echo "=========================================="
echo "Callers: ${WGS_CALLERS}"
echo "Coverages: ${COVERAGES_WGS[*]}"
echo ""

for COV in "${COVERAGES_WGS[@]}"; do
    echo "--- Coverage: ${COV}x ---"

    # Temporary directory for per-caller variant lists
    TMP_DIR=$(mktemp -d)
    trap "rm -rf ${TMP_DIR}" EXIT

    CALLER_LIST=()
    ALIAS_LIST=()

    for caller in ${WGS_CALLERS}; do
        alias=$(caller_alias "${caller}")
        CALLER_DIR="${VCFEVAL_ROOT}/${COV}x/${caller}"

        if [[ ! -d "${CALLER_DIR}" ]]; then
            echo "  SKIP: ${alias} — no vcfeval output"
            continue
        fi

        CALLER_LIST+=("${caller}")
        ALIAS_LIST+=("${alias}")

        # Extract TP variants (called=1, truth=1)
        echo "  Extracting TP for ${alias}..."
        extract_variant_ids "${CALLER_DIR}/tp.vcf.gz" | sort -u > "${TMP_DIR}/${alias}_tp.txt"

        # Extract FP variants (called=1, truth=0)
        echo "  Extracting FP for ${alias}..."
        extract_variant_ids "${CALLER_DIR}/fp.vcf.gz" | sort -u > "${TMP_DIR}/${alias}_fp.txt"

        # Extract FN variants from tp-baseline (called=0, truth=1)
        echo "  Extracting FN for ${alias}..."
        extract_variant_ids "${CALLER_DIR}/fn.vcf.gz" | sort -u > "${TMP_DIR}/${alias}_fn.txt"

        # All variants this caller reported (TP + FP)
        cat "${TMP_DIR}/${alias}_tp.txt" "${TMP_DIR}/${alias}_fp.txt" | sort -u > "${TMP_DIR}/${alias}_called.txt"
    done

    if [[ ${#CALLER_LIST[@]} -lt 2 ]]; then
        echo "  SKIP: fewer than 2 callers found"
        continue
    fi

    # Collect all unique variant IDs across all callers
    echo "  Merging variant IDs..."

    # All truth-positive variants (union of TP + FN across callers)
    cat "${TMP_DIR}"/*_tp.txt "${TMP_DIR}"/*_fn.txt | sort -u > "${TMP_DIR}/all_truth.txt"

    # All called variants (union of TP + FP across callers)
    cat "${TMP_DIR}"/*_called.txt | sort -u > "${TMP_DIR}/all_called.txt"

    # All unique variant IDs
    cat "${TMP_DIR}/all_truth.txt" "${TMP_DIR}/all_called.txt" | sort -u > "${TMP_DIR}/all_variants.txt"

    N_VARIANTS=$(wc -l < "${TMP_DIR}/all_variants.txt")
    echo "  Total unique variants: ${N_VARIANTS}"

    # =========================================================================
    # 1) Build concordance matrix TSV
    # =========================================================================
    CONCORDANCE_TSV="${OUT_DIR}/variant_concordance_${COV}x.tsv"
    echo "  Building concordance matrix..."

    # Header
    printf 'variant_id\ttruth'  > "${CONCORDANCE_TSV}"
    for alias in "${ALIAS_LIST[@]}"; do
        printf '\t%s' "${alias}" >> "${CONCORDANCE_TSV}"
    done
    printf '\n' >> "${CONCORDANCE_TSV}"

    # For each variant, determine truth and per-caller call status
    while IFS= read -r var_id; do
        # Is it in truth set?
        truth=0
        if grep -qFx "${var_id}" "${TMP_DIR}/all_truth.txt"; then
            truth=1
        fi

        printf '%s\t%d' "${var_id}" "${truth}" >> "${CONCORDANCE_TSV}"

        for alias in "${ALIAS_LIST[@]}"; do
            called=0
            if grep -qFx "${var_id}" "${TMP_DIR}/${alias}_called.txt"; then
                called=1
            fi
            printf '\t%d' "${called}" >> "${CONCORDANCE_TSV}"
        done
        printf '\n' >> "${CONCORDANCE_TSV}"
    done < "${TMP_DIR}/all_variants.txt"

    echo "  Written: ${CONCORDANCE_TSV}"

    # =========================================================================
    # 2) Upset plot long-format: FP (callset) — variants called but not in truth
    # =========================================================================
    UPSET_CALLSET="${OUT_DIR}/upset_long_callset_${COV}x.csv"
    echo "  Building upset FP (callset) data..."
    > "${UPSET_CALLSET}"

    for alias in "${ALIAS_LIST[@]}"; do
        while IFS= read -r var_id; do
            echo "${alias},${var_id}" >> "${UPSET_CALLSET}"
        done < "${TMP_DIR}/${alias}_fp.txt"
    done

    echo "  Written: ${UPSET_CALLSET}"

    # =========================================================================
    # 3) Upset plot long-format: FN (baseline) — truth variants missed
    # =========================================================================
    UPSET_BASELINE="${OUT_DIR}/upset_long_baseline_${COV}x.csv"
    echo "  Building upset FN (baseline) data..."
    > "${UPSET_BASELINE}"

    for alias in "${ALIAS_LIST[@]}"; do
        while IFS= read -r var_id; do
            echo "${alias},${var_id}" >> "${UPSET_BASELINE}"
        done < "${TMP_DIR}/${alias}_fn.txt"
    done

    echo "  Written: ${UPSET_BASELINE}"

    # Cleanup temp
    rm -rf "${TMP_DIR}"
    trap - EXIT

    echo ""
done

echo "=========================================="
echo " Variant concordance tables complete"
echo "=========================================="
echo "Files per coverage:"
echo "  variant_concordance_{cov}x.tsv  — binary matrix"
echo "  upset_long_callset_{cov}x.csv   — FP long format"
echo "  upset_long_baseline_{cov}x.csv  — FN long format"
