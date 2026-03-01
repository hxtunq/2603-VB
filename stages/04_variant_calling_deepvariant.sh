#!/bin/bash
#===============================================================================
# STEP 04: Variant Calling - DeepVariant 1.8.0 (via Docker)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/../scripts/helper_functions.sh"

CALLER="deepvariant"
log_info "===== STEP 04: ${CALLER} (v${DEEPVARIANT_VERSION}, Docker) ====="

# Input
source "${PREPROC_DIR}/bam_path.sh"

OUT_DIR="${VARIANT_DIR}/${CALLER}"
ensure_dir "${OUT_DIR}/intermediate"

#-------------------------------------------------------------------------------
# 1. Prepare paths for Docker
#-------------------------------------------------------------------------------
# Docker needs absolute paths
ABS_REF_DIR=$(cd "${REF_DIR}" && pwd)
ABS_PREPROC_DIR=$(cd "${PREPROC_DIR}" && pwd)
ABS_OUT_DIR=$(cd "${OUT_DIR}" && pwd)

BAM_BASENAME=$(basename "${FINAL_BAM}")
REF_BASENAME=$(basename "${REF_FASTA}")

#-------------------------------------------------------------------------------
# 2. Run DeepVariant via Docker
#-------------------------------------------------------------------------------
log_info "Running DeepVariant via Docker..."
log_info "  Image: ${DEEPVARIANT_IMAGE}"
log_info "  Model: WGS"
start_timer

DOCKER_USER="$(id -u)":"$(id -g)"

docker run \
    --rm \
    --user "${DOCKER_USER}" \
    --cpus "${THREADS}" \
    --memory "${MAX_MEMORY}" \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_PREPROC_DIR}:/input:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    ${DEEPVARIANT_IMAGE} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/ref/${REF_BASENAME}" \
    --reads="/input/${BAM_BASENAME}" \
    --sample_name="${SAMPLE_NAME}" \
    --output_vcf="/output/${PREFIX}_${CALLER}_raw.vcf.gz" \
    --output_gvcf="/output/${PREFIX}_${CALLER}.g.vcf.gz" \
    --intermediate_results_dir="/output/intermediate" \
    --num_shards="${THREADS}" \
    2>&1 | tee "${LOG_DIR}/${CALLER}.log"

end_timer "04_${CALLER}"
log_info "DeepVariant completed"

#-------------------------------------------------------------------------------
# 3. Process output
#-------------------------------------------------------------------------------
log_info "Processing DeepVariant output..."

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf.gz"
PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"

# Ensure sorted/bgzip for tabix (DeepVariant output can be unsorted)
SORT_TMP="${OUT_DIR}/intermediate/bcftools_sort"
ensure_dir "${SORT_TMP}"
SORTED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.sorted.vcf.gz"
bcftools sort "${RAW_VCF}" -Oz -o "${SORTED_VCF}" -T "${SORT_TMP}"
mv -f "${SORTED_VCF}" "${RAW_VCF}"

# Index
tabix -f -p vcf "${RAW_VCF}"

# DeepVariant outputs mostly PASS, but filter to be safe
bcftools view -f "PASS,." "${RAW_VCF}" -Oz -o "${PASS_VCF}"
tabix -f -p vcf "${PASS_VCF}"

# Split by type
bcftools view -v snps "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 4. Normalize for benchmarking
#-------------------------------------------------------------------------------
log_info "Normalizing variants for benchmarking..."

NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
"${SCRIPT_DIR}/../scripts/normalize_vcf.sh" "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
ensure_dir "$(dirname "${TRUTH_NORM}")"
"${SCRIPT_DIR}/../scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"

#-------------------------------------------------------------------------------
# 5. Stats
#-------------------------------------------------------------------------------
bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${PASS_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${PASS_VCF}" | wc -l)

log_info "Results: $((N_SNP + N_INDEL)) variants (${N_SNP} SNPs, ${N_INDEL} INDELs)"

# Cleanup intermediate files
rm -rf "${OUT_DIR}/intermediate"

log_info "===== ${CALLER} Complete ====="
