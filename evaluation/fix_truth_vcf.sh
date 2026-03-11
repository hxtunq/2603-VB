#!/usr/bin/env bash
# fix_truth_vcf.sh — Add sample column (GT=1/1) to simuG sites-only truth VCF.
#
# simuG outputs sites-only VCFs (no FORMAT/SAMPLE columns).
# Because the simulation is haploid (ART reads come from one mutated genome),
# all variant positions have 100% ALT reads → callers report GT=1/1.
# This script adds FORMAT=GT and SAMPLE=1/1 so hap.py can match genotypes.
#
# Usage:
#   bash evaluation/fix_truth_vcf.sh
#
# Prerequisite: Step 3.1 of workflow.md must have already produced the
# sites-only truth VCF at ${TRUTH_VCF}.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

SITES_VCF="${TRUTH_VCF}"     # e.g. data/simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz
FIXED_VCF="${SITES_VCF%.vcf.gz}_fixed.vcf.gz"

if [[ ! -f "${SITES_VCF}" ]]; then
    echo "ERROR: Truth VCF not found: ${SITES_VCF}"
    exit 1
fi

echo "=== Fixing truth VCF: adding GT=1/1 sample column ==="
echo "  Input:  ${SITES_VCF}"
echo "  Output: ${FIXED_VCF}"

# Decompress → add FORMAT + SAMPLE columns → compress + index
zcat "${SITES_VCF}" | awk '
BEGIN { OFS="\t" }
/^##/ { print; next }
/^#CHROM/ {
    # Add FORMAT header before #CHROM line
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    # Append FORMAT and SAMPLE columns to header
    print $0, "FORMAT", "SIMULATED_SAMPLE"
    next
}
{
    # Data lines: append GT=1/1 (haploid sim → 100% ALT reads → hom ALT)
    print $0, "GT", "1/1"
}
' | bgzip -c > "${FIXED_VCF}"

tabix -f -p vcf "${FIXED_VCF}"

# Verify
echo ""
echo "=== Verification ==="
echo "Header sample check:"
bcftools query -l "${FIXED_VCF}"
echo ""
echo "First 5 variants (GT column):"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "${FIXED_VCF}" | head -5 || true
echo ""

NSNP=$(bcftools view -v snps "${FIXED_VCF}" | grep -cv '^#' || true)
NINDEL=$(bcftools view -v indels "${FIXED_VCF}" | grep -cv '^#' || true)
echo "SNP count:   ${NSNP}"
echo "INDEL count: ${NINDEL}"

# Replace original with fixed version
echo ""
echo "Replacing original truth VCF..."
mv "${FIXED_VCF}" "${SITES_VCF}"
mv "${FIXED_VCF}.tbi" "${SITES_VCF}.tbi"
echo "Done: ${SITES_VCF} now has GT=1/1 sample column."
