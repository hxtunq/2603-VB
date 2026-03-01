#!/bin/bash
#===============================================================================
# Helper script: Rename chromosomes from "22" to "chr22" format
# Usage: ./rename_chromosomes.sh input.vcf.gz output.vcf.gz
#===============================================================================

set -euo pipefail

INPUT_VCF="$1"
OUTPUT_VCF="$2"

# Create chromosome mapping file
CHR_MAP=$(mktemp)
for i in {1..22} X Y M MT; do
    echo -e "${i}\tchr${i}" >> "${CHR_MAP}"
done
# Handle MT -> chrM
echo -e "MT\tchrM" >> "${CHR_MAP}"

bcftools annotate --rename-chrs "${CHR_MAP}" "${INPUT_VCF}" -Oz -o "${OUTPUT_VCF}"

# Index output
tabix -f -p vcf "${OUTPUT_VCF}"

# Cleanup
rm -f "${CHR_MAP}"

echo "Done: ${OUTPUT_VCF}"
