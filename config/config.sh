#!/bin/bash
#===============================================================================
# CONFIG: Cấu hình chung cho pipeline variant calling benchmarking
# Theo GATK Best Practices (nf-core/sarek) — Multi-coverage, WGS + WES
#===============================================================================

#-------------------------------------------------------------------------------
# SYSTEM RESOURCES (auto-detect with safe caps)
#-------------------------------------------------------------------------------
CPU_DETECTED=$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 4)
if [[ -z "${CPU_DETECTED}" || "${CPU_DETECTED}" -lt 1 ]]; then
    CPU_DETECTED=4
fi

# Keep one core free for system responsiveness, and cap at 8 for local benchmark fairness.
THREADS=$((CPU_DETECTED - 1))
if [[ "${THREADS}" -lt 1 ]]; then
    THREADS=1
elif [[ "${THREADS}" -gt 8 ]]; then
    THREADS=8
fi
export THREADS

# Cap memory usage to machine limit (14G requested by project constraints).
export MAX_MEMORY="14G"
export JAVA_HEAP="12G"
export JAVA_OPTS="-Xmx${JAVA_HEAP} -XX:ParallelGCThreads=2"

#-------------------------------------------------------------------------------
# DIRECTORY PATHS
#-------------------------------------------------------------------------------
CONFIG_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PROJECT_DIR="${PROJECT_DIR:-$(cd "${CONFIG_DIR}/.." && pwd)}"
export DATA_DIR="${PROJECT_DIR}/data"
export REF_DIR="${DATA_DIR}/reference"
export SIM_DIR="${DATA_DIR}/simulated"
export RESULTS_DIR="${PROJECT_DIR}/results"
export LOG_DIR="${PROJECT_DIR}/logs"

export PREPROC_DIR="${RESULTS_DIR}/preprocessing"
export VARIANT_DIR="${RESULTS_DIR}/variants"
export BENCH_DIR="${RESULTS_DIR}/benchmarks"
export EVAL_DIR="${RESULTS_DIR}/eval"
export FIGURE_DIR="${RESULTS_DIR}/plots"
export METRICS_DIR="${RESULTS_DIR}/final_metrics"

#-------------------------------------------------------------------------------
# REFERENCE GENOME
#-------------------------------------------------------------------------------
export GENOME_VERSION="hg38"
export CHR_TO_USE="chr22"
export REF_FASTA="${REF_DIR}/${CHR_TO_USE}.fa"
export REF_DICT="${REF_DIR}/${CHR_TO_USE}.dict"
export REF_FAI="${REF_DIR}/${CHR_TO_USE}.fa.fai"

# Known sites for BQSR (GATK Bundle)
export DBSNP="${REF_DIR}/dbsnp138.hg38.chr22.vcf.gz"
export KNOWN_INDELS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz"
export KNOWN_SNPS="${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz"

#-------------------------------------------------------------------------------
# SIMULATION PARAMETERS (simuG + ART)
#-------------------------------------------------------------------------------
# simuG — realistic for chr22 (~50 Mb, ~20K true variants)
# NOTE: In simuG random SNP mode (-snp_count), set -titv_ratio explicitly in
# the calling command; omitting it falls back to the tool default of 0.5.
export SNP_COUNT=18000          # ~18K SNPs (realistic for chr22)
export INDEL_COUNT=2500         # ~2.5K indels (SNP:INDEL ~ 7:1)
export INDEL_MIN_LEN=1          # Minimum indel length
export INDEL_MAX_LEN=30         # Up to 30bp to test longer indels

# ART Illumina — multi-coverage design
# WGS coverages: low to medium (like GIAB WGS 22-37x)
# Higher coverages: for WES-like evaluation (GIAB WES 183-249x)
export COVERAGES_WGS=(10 20 30 50)
export COVERAGES_WES=(50 100 200)
export COVERAGES_ALL=(10 20 30 50 100 200)

export READ_LENGTH=150
export FRAGMENT_MEAN=350
export FRAGMENT_SD=50
export ART_PLATFORM="HS25"
export SEED=42

#-------------------------------------------------------------------------------
# QUALITY THRESHOLDS
#-------------------------------------------------------------------------------
export MIN_BASE_QUALITY=20
export MIN_MAPPING_QUALITY=20
export MIN_READ_LENGTH=50

#-------------------------------------------------------------------------------
# SAMPLE INFO
#-------------------------------------------------------------------------------
export SAMPLE_NAME="SIMULATED_SAMPLE"
export PREFIX="${SAMPLE_NAME}_${CHR_TO_USE}"
export READ_GROUP="@RG\\tID:${SAMPLE_NAME}\\tSM:${SAMPLE_NAME}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"

#-------------------------------------------------------------------------------
# BED FILES — Exome targets + stratification
#-------------------------------------------------------------------------------
# Exome target BED (Agilent SureSelect V6, extracted for chr22)
export EXOME_BED="${REF_DIR}/Exome-Agilent_V6.chr22.bed"

# CDS-only BED (canonical protein-coding, extracted for chr22)
export CDS_BED="${REF_DIR}/CDS-canonical.chr22.bed"

# Mappability tricky regions (low-mappability, extracted for chr22)
export MAPPABILITY_BED="${REF_DIR}/umap_k100_mappability.chr22.bed.gz"

# Full-chromosome high-confidence region (for simulated data = non-N regions)
export TRUTH_BED="${REF_DIR}/${CHR_TO_USE}_highconf.bed"

# Stratification TSV for hap.py (list of BED files)
export STRATIFICATION_TSV="${REF_DIR}/stratification_chr22.tsv"

#-------------------------------------------------------------------------------
# TOOL VERSIONS
#-------------------------------------------------------------------------------
# GATK 4.6.2.0 (system-installed)
# FreeBayes 1.3.10 (system-installed)

#-------------------------------------------------------------------------------
# DOCKER IMAGES
#-------------------------------------------------------------------------------
export DEEPVARIANT_VERSION="1.9.0"
export DEEPVARIANT_IMAGE="google/deepvariant:${DEEPVARIANT_VERSION}"
export STRELKA2_IMAGE="quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
export MANTA_IMAGE="quay.io/biocontainers/manta:1.6.0--h9ee0642_2"
export RTG_IMAGE="biocontainers/rtg-tools:3.12.1--hdfd78af_0"
export HAPPY_IMAGE="jmcdani20/hap.py:v0.3.12"

# Sentieon DNAscope (optional — requires SENTIEON_LICENSE)
export SENTIEON_VERSION="${SENTIEON_VERSION:-202503.02}"

# === Sentieon DNAscope Model Bundles ===
# Download from: https://github.com/Sentieon/sentieon-models
# sentieon-cli dnascope expects -m/--model_bundle to be the .bundle archive file.
# Do not point these variables at unpacked directories.
# Override DNASCOPE_WGS_MODEL / DNASCOPE_WES_MODEL if you download newer bundles.
# WGS: DNAscopeIlluminaWGS2.0.bundle
# WES: DNAscopeIlluminaWES2.0.bundle
export DNASCOPE_WGS_MODEL="${DNASCOPE_WGS_MODEL:-${REF_DIR}/models/DNAscopeIlluminaWGS2.0.bundle}"
export DNASCOPE_WES_MODEL="${DNASCOPE_WES_MODEL:-${REF_DIR}/models/DNAscopeIlluminaWES2.0.bundle}"

# PCR-free library prep (affects DNAscope indel model)
export PCRFREE="${PCRFREE:-true}"

#-------------------------------------------------------------------------------
# BENCHMARKING
#-------------------------------------------------------------------------------
export RTG_SDF="${REF_DIR}/${CHR_TO_USE}.sdf"

#-------------------------------------------------------------------------------
# VARIANT CALLER PARAMETERS
#-------------------------------------------------------------------------------
# GATK HaplotypeCaller (default GATK 4.x = 10)
export GATK_STAND_CALL_CONF=10

# FreeBayes
export FB_MIN_ALT_COUNT=3
export FB_MIN_ALT_FRACTION=0.2

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------
export TRUTH_VCF="${SIM_DIR}/${PREFIX}_truth.vcf.gz"

echo "[CONFIG] Loaded successfully"
echo "[CONFIG] Project: ${PROJECT_DIR}"
echo "[CONFIG] Chromosome: ${CHR_TO_USE}"
echo "[CONFIG] WGS coverages: ${COVERAGES_WGS[*]}"
echo "[CONFIG] WES coverages: ${COVERAGES_WES[*]}"
