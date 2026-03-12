#!/usr/bin/env bash
# =============================================================================
# Download GIAB v3.3 stratification BED files for repeat/mappability/segdup.
#
# Downloads from GIAB FTP, filters for chr22, and saves to ${REF_DIR}/repeats/.
# Appends new entries to stratification_chr22.tsv.
#
# BED files downloaded:
#   1. Tandem repeats + homopolymers (LowComplexity)
#   2. Low mappability (all)
#   3. Segmental duplications (>10kb)
#
# Usage:
#   bash evaluation/prepare_giab_strata.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

REPEATS_DIR="${REPEATS_DIR:-${REF_DIR}/repeats}"
mkdir -p "${REPEATS_DIR}"

GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38"

echo "=========================================="
echo " Downloading GIAB v3.3 Stratification BEDs"
echo "=========================================="

# -----------------------------------------------------------------------------
# Helper: download, gunzip, filter chr22, sort, merge
# -----------------------------------------------------------------------------
download_and_filter() {
    local URL="$1"
    local OUT="$2"
    local LABEL="$3"

    if [[ -f "${OUT}" ]] && [[ -s "${OUT}" ]]; then
        echo "  Already exists: ${OUT} ($(wc -l < "${OUT}") intervals)"
        return 0
    fi

    local TMP="${OUT}.tmp.gz"
    echo "  Downloading ${LABEL}..."
    wget -q -O "${TMP}" "${URL}" || { echo "  ERROR: Download failed for ${LABEL}"; rm -f "${TMP}"; return 1; }

    echo "  Filtering for ${CHR_TO_USE}..."
    zcat "${TMP}" \
        | awk -v chr="${CHR_TO_USE}" '$1 == chr' \
        | sort -k1,1 -k2,2n \
        > "${OUT}"
    rm -f "${TMP}"

    local N_INTERVALS
    N_INTERVALS=$(wc -l < "${OUT}")
    if [[ "${N_INTERVALS}" -eq 0 ]]; then
        echo "  WARNING: No intervals found for ${CHR_TO_USE} in ${LABEL}"
    else
        echo "  Created: ${OUT} (${N_INTERVALS} intervals)"
    fi
}

# =============================================================================
# 1) Tandem Repeats + Homopolymers
# =============================================================================

echo ""
echo "--- 1. Tandem Repeats + Homopolymers ---"
download_and_filter \
    "${GIAB_BASE}/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz" \
    "${REPEATS_DIR}/tandem_repeat_homopolymer.bed" \
    "AllTandemRepeatsandHomopolymers"

# =============================================================================
# 2) Low Mappability (all)
# =============================================================================

echo ""
echo "--- 2. Low Mappability ---"
download_and_filter \
    "${GIAB_BASE}/mappability/GRCh38_lowmappabilityall.bed.gz" \
    "${REPEATS_DIR}/low_mappability.bed" \
    "LowMappabilityAll"

# =============================================================================
# 3) Segmental Duplications (>10kb)
# =============================================================================

echo ""
echo "--- 3. Segmental Duplications ---"
download_and_filter \
    "${GIAB_BASE}/SegmentalDuplications/GRCh38_segdups_gt10kb.bed.gz" \
    "${REPEATS_DIR}/segdup.bed" \
    "SegDup_gt10kb"

# =============================================================================
# 4) Append to stratification TSV (idempotent)
# =============================================================================

echo ""
echo "--- 4. Updating stratification TSV ---"

STRAT_TSV="${REF_DIR}/stratification_chr22.tsv"

# Create if not exists
touch "${STRAT_TSV}"

append_if_missing() {
    local NAME="$1"
    local REL_PATH="$2"
    if ! grep -q "^${NAME}	" "${STRAT_TSV}" 2>/dev/null; then
        echo "${NAME}	${REL_PATH}" >> "${STRAT_TSV}"
        echo "  Added: ${NAME}"
    else
        echo "  Already in TSV: ${NAME}"
    fi
}

append_if_missing "tandem_repeat_homopolymer" "repeats/tandem_repeat_homopolymer.bed"
append_if_missing "low_mappability"           "repeats/low_mappability.bed"
append_if_missing "segdup_gt10kb"             "repeats/segdup.bed"

echo ""
echo "=========================================="
echo " GIAB stratification download complete"
echo "=========================================="
echo "  BED files: ${REPEATS_DIR}/"
echo "  TSV:       ${STRAT_TSV} ($(wc -l < "${STRAT_TSV}") entries)"
