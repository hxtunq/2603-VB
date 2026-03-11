#!/usr/bin/env bash
# =============================================================================
# Download and filter ClinVar data for hg38 / chr22 benchmarking.
#
# 1. Downloads latest ClinVar VCF (hg38) from NCBI FTP
# 2. Filters to chr22
# 3. Keeps Pathogenic + Likely_pathogenic (non-conflicting)
# 4. Creates BED file for stratification
# 5. Appends to stratification_chr22.tsv
#
# Usage:
#   bash evaluation/prepare_clinvar.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_RAW="${REF_DIR}/clinvar_hg38.vcf.gz"
CLINVAR_CHR22="${REF_DIR}/clinvar_chr22.vcf.gz"
CLINVAR_FILTERED="${REF_DIR}/clinvar_chr22_pathogenic.vcf.gz"
CLINVAR_BED="${REF_DIR}/stratification/clinvar_chr22_pathogenic.bed"
STRAT_TSV="${REF_DIR}/stratification_chr22.tsv"

echo "=========================================="
echo " ClinVar Data Preparation (hg38 / chr22)"
echo "=========================================="

# ---- Pre-flight ----
command -v bcftools &>/dev/null || { echo "ERROR: 'bcftools' not found in PATH"; exit 1; }
command -v bgzip &>/dev/null   || { echo "ERROR: 'bgzip' not found in PATH"; exit 1; }
command -v tabix &>/dev/null   || { echo "ERROR: 'tabix' not found in PATH"; exit 1; }

# =============================================================================
# 1) Download ClinVar VCF
# =============================================================================

echo ""
echo "--- Step 1: Download ClinVar VCF ---"

if [[ ! -f "${CLINVAR_RAW}" ]]; then
    echo "  Downloading ClinVar (hg38)..."
    wget -q -O "${CLINVAR_RAW}" "${CLINVAR_URL}"
    wget -q -O "${CLINVAR_RAW}.tbi" "${CLINVAR_URL}.tbi"
    echo "  Downloaded: ${CLINVAR_RAW}"
else
    echo "  Already exists: ${CLINVAR_RAW}"
fi

# =============================================================================
# 2) Filter to chr22
# =============================================================================

echo ""
echo "--- Step 2: Filter to chr22 ---"

if [[ ! -f "${CLINVAR_CHR22}" ]]; then
    echo "  Extracting chr22..."
    # ClinVar uses "chr22" for GRCh38 since 2023; older versions use "22"
    # Try chr22 first, fall back to 22
    CHR_FOUND=$(bcftools view -h "${CLINVAR_RAW}" | grep "^##contig" | grep -oP 'ID=(chr22|22)' | head -1 || echo "")

    if [[ "${CHR_FOUND}" == "ID=chr22" ]]; then
        bcftools view -r chr22 "${CLINVAR_RAW}" | bgzip -c > "${CLINVAR_CHR22}"
    elif [[ "${CHR_FOUND}" == "ID=22" ]]; then
        echo "  ClinVar uses numeric chromosomes — extracting '22' and renaming to 'chr22'..."
        bcftools view -r 22 "${CLINVAR_RAW}" \
            | bcftools annotate --rename-chrs <(echo "22 chr22") \
            | bgzip -c > "${CLINVAR_CHR22}"
    else
        echo "  WARNING: Could not detect chr22 format; trying both..."
        bcftools view -r chr22,22 "${CLINVAR_RAW}" \
            | sed 's/^22\t/chr22\t/' \
            | bgzip -c > "${CLINVAR_CHR22}"
    fi
    tabix -f -p vcf "${CLINVAR_CHR22}"

    N_VARIANTS=$(bcftools view -H "${CLINVAR_CHR22}" | wc -l)
    echo "  chr22 variants: ${N_VARIANTS}"
else
    echo "  Already exists: ${CLINVAR_CHR22}"
fi

# =============================================================================
# 3) Filter to Pathogenic / Likely_pathogenic (non-conflicting)
# =============================================================================

echo ""
echo "--- Step 3: Filter pathogenic variants ---"

if [[ ! -f "${CLINVAR_FILTERED}" ]]; then
    echo "  Filtering Pathogenic + Likely_pathogenic (non-conflicting)..."

    # ClinVar INFO field: CLNSIG contains clinical significance
    # Filter logic:
    #   - Keep: CLNSIG contains "Pathogenic" or "Likely_pathogenic"
    #   - Exclude: CLNSIG contains "Conflicting"
    #   - Also exclude CLNSIGCONF (conflicting interpretations field, newer ClinVar)
    bcftools view -H "${CLINVAR_CHR22}" \
        | awk -F'\t' '
        {
            info = $8
            # Extract CLNSIG value
            clnsig = ""
            if (match(info, /CLNSIG=[^;]+/)) {
                clnsig = substr(info, RSTART+7, RLENGTH-7)
            }

            # Skip if conflicting
            if (clnsig ~ /Conflicting/) next
            if (info ~ /CLNSIGCONF=/) next

            # Keep if Pathogenic or Likely_pathogenic
            if (clnsig ~ /Pathogenic/ || clnsig ~ /Likely_pathogenic/) print $0
        }' \
        | cat <(bcftools view -h "${CLINVAR_CHR22}") - \
        | bgzip -c > "${CLINVAR_FILTERED}"
    tabix -f -p vcf "${CLINVAR_FILTERED}"

    N_PATH=$(bcftools view -H "${CLINVAR_FILTERED}" | wc -l)
    echo "  Pathogenic/Likely_pathogenic variants on chr22: ${N_PATH}"

    # Count by type
    N_SNP=$(bcftools view -H "${CLINVAR_FILTERED}" | awk 'length($4)==1 && length($5)==1' | wc -l)
    N_INDEL=$(bcftools view -H "${CLINVAR_FILTERED}" | awk 'length($4)!=length($5) || length($4)>1' | wc -l)
    echo "    SNPs:   ${N_SNP}"
    echo "    INDELs: ${N_INDEL}"
else
    echo "  Already exists: ${CLINVAR_FILTERED}"
fi

# =============================================================================
# 4) Create BED file for stratification
# =============================================================================

echo ""
echo "--- Step 4: Creating ClinVar BED ---"

mkdir -p "$(dirname "${CLINVAR_BED}")"

if [[ ! -f "${CLINVAR_BED}" ]]; then
    bcftools query -f '%CHROM\t%POS0\t%END\n' "${CLINVAR_FILTERED}" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "${CLINVAR_BED}"
    echo "  Created: ${CLINVAR_BED} ($(wc -l < "${CLINVAR_BED}") intervals)"
else
    echo "  Already exists: ${CLINVAR_BED}"
fi

# =============================================================================
# 5) Append to stratification TSV
# =============================================================================

echo ""
echo "--- Step 5: Updating stratification TSV ---"

if [[ -f "${STRAT_TSV}" ]]; then
    # Remove old ClinVar entry if present
    grep -v "clinvar" "${STRAT_TSV}" > "${STRAT_TSV}.tmp" || true
    mv "${STRAT_TSV}.tmp" "${STRAT_TSV}"
fi

echo "clinvar_pathogenic	stratification/clinvar_chr22_pathogenic.bed" >> "${STRAT_TSV}"
echo "  Updated: ${STRAT_TSV}"

echo ""
echo "=========================================="
echo " ClinVar preparation complete"
echo "=========================================="
echo "  Filtered VCF: ${CLINVAR_FILTERED}"
echo "  BED file:     ${CLINVAR_BED}"
echo ""
echo "Next: bash evaluation/eval_happy.sh (with stratification)"
