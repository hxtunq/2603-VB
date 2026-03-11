#!/usr/bin/env bash
# fix_truth_vcf.sh - Ensure the simuG truth VCF has exactly one GT sample column.
#
# simuG emits sites-only VCFs by default. This script adds FORMAT=GT and one
# sample column with GT=1/1. If the input already has sample columns, the script
# rewrites the VCF to keep only the first sample column so reruns stay safe.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

SITES_VCF="${TRUTH_VCF}"     # e.g. data/simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz
FIXED_VCF="${SITES_VCF%.vcf.gz}_fixed.vcf.gz"

if [[ ! -f "${SITES_VCF}" ]]; then
    echo "ERROR: Truth VCF not found: ${SITES_VCF}"
    exit 1
fi

echo "=== Fixing truth VCF: ensuring one GT sample column ==="
echo "  Input:  ${SITES_VCF}"
echo "  Output: ${FIXED_VCF}"

sample_count=$(bcftools query -l "${SITES_VCF}" | awk 'NF {c++} END {print c+0}')

if (( sample_count > 0 )); then
    echo "  Input already has ${sample_count} sample(s); keeping only the first sample column."

    zcat "${SITES_VCF}" | awk -v sample="${SAMPLE_NAME}" '
    BEGIN { OFS="\t"; ff=0; gt=0 }
    /^##fileformat=/ { ff=1; print; next }
    /^##FORMAT=<ID=GT,/ { gt=1; print; next }
    /^##/ { print; next }
    /^#CHROM/ {
        if (ff == 0) {
            print "##fileformat=VCFv4.2"
        }
        if (gt == 0) {
            print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        }
        print $1, $2, $3, $4, $5, $6, $7, $8, "FORMAT", sample
        next
    }
    {
        if (NF < 10) {
            printf "ERROR: expected FORMAT + sample columns at %s:%s but found %d fields\n", $1, $2, NF > "/dev/stderr"
            exit 1
        }
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
    }
    ' | bgzip -c > "${FIXED_VCF}"
else
    echo "  Input is sites-only; adding GT=1/1 sample column."

    zcat "${SITES_VCF}" | awk -v sample="${SAMPLE_NAME}" '
    BEGIN { OFS="\t"; ff=0 }
    /^##fileformat=/ { ff=1; print; next }
    /^##/ { print; next }
    /^#CHROM/ {
        if (ff == 0) {
            print "##fileformat=VCFv4.2"
        }
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print $0, "FORMAT", sample
        next
    }
    {
        print $0, "GT", "1/1"
    }
    ' | bgzip -c > "${FIXED_VCF}"
fi

tabix -f -p vcf "${FIXED_VCF}"

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

echo ""
echo "Replacing original truth VCF..."
mv "${FIXED_VCF}" "${SITES_VCF}"
mv "${FIXED_VCF}.tbi" "${SITES_VCF}.tbi"
echo "Done: ${SITES_VCF} now has one GT sample column."
