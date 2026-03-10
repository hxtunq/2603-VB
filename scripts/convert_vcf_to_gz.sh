#!/bin/bash
# Convert all plain .vcf output files to .vcf.gz (bgzip + tabix)
# Usage: bash scripts/convert_vcf_to_gz.sh
# Run from: variant-calling-benchmark/

set -euo pipefail

source "$(dirname "$0")/../config/config.sh"

CONVERTED=0
SKIPPED=0

echo "=== Converting .vcf to .vcf.gz in ${VARIANT_DIR} ==="

# Find all final output VCFs (SS_chr22_*.vcf) but exclude intermediate files
find "${VARIANT_DIR}" -name "SS_${CHR_TO_USE}_*.vcf" -type f | sort | while read -r vcf; do
    gz="${vcf}.gz"

    # Skip if .vcf.gz already exists
    if [[ -f "${gz}" ]]; then
        echo "  SKIP (gz exists): $(basename "${vcf}")"
        ((SKIPPED++)) || true
        continue
    fi

    echo "  CONVERT: ${vcf}"
    bgzip -c "${vcf}" > "${gz}"
    tabix -f -p vcf "${gz}"
    ((CONVERTED++)) || true
done

echo ""
echo "Done. Converted: ${CONVERTED}, Skipped: ${SKIPPED}"
echo "Original .vcf files are preserved. Remove them manually if no longer needed:"
echo "  find ${VARIANT_DIR} -name 'SS_${CHR_TO_USE}_*.vcf' -type f -delete"
