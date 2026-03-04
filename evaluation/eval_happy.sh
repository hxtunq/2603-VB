#!/usr/bin/env bash
# Evaluate variant calls with hap.py — Multi-coverage, WGS + WES
# Based on: caller_benchmark-main/eval_happy.sh
# Adapted for multi-coverage simulated benchmark with stratification

set -e

source "$(dirname "$0")/../config/config.sh"

# ============================================================
# Caller definitions — maps to VCF file patterns
# ============================================================
WGS_CALLERS="gatk deepvariant strelka2 freebayes dnascope"
WES_CALLERS="gatk_wes deepvariant_wes strelka2_wes freebayes_wes dnascope_wes"

get_vcf() {
    local caller=$1 variant_dir=$2
    case $caller in
        gatk)              echo "${variant_dir}/gatk/${PREFIX}_HC_HARDFILTER.vcf" ;;
        deepvariant)       echo "${variant_dir}/deepvariant/${PREFIX}_DV_STANDART.vcf" ;;
        strelka2)          echo "${variant_dir}/strelka2/${PREFIX}_STRELKA_STANDART.vcf.gz" ;;
        freebayes)         echo "${variant_dir}/freebayes/${PREFIX}_FB_STANDART.vcf" ;;
        dnascope)          echo "${variant_dir}/dnascope/${PREFIX}_DNASCOPE.vcf" ;;
        gatk_wes)          echo "${variant_dir}/gatk_wes/${PREFIX}_HC_HARDFILTER.vcf" ;;
        deepvariant_wes)   echo "${variant_dir}/deepvariant_wes/${PREFIX}_DV_STANDART.vcf" ;;
        strelka2_wes)      echo "${variant_dir}/strelka2_wes/${PREFIX}_STRELKA_STANDART.vcf.gz" ;;
        freebayes_wes)     echo "${variant_dir}/freebayes_wes/${PREFIX}_FB_STANDART.vcf" ;;
        dnascope_wes)      echo "${variant_dir}/dnascope_wes/${PREFIX}_DNASCOPE.vcf" ;;
    esac
}

# ============================================================
# Evaluation function
# ============================================================
function make_comparison()
{
    local CURR_VCF="$1"
    local CALLER_NAME="$2"
    local COV="$3"
    local MODE="$4"  # WGS or WES

    local EVAL_OUT="${EVAL_DIR}/${COV}x/${MODE}/${CALLER_NAME}_eval_data"
    mkdir -p "${EVAL_OUT}"

    echo "Evaluating: ${CALLER_NAME} @ ${COV}x (${MODE})"
    echo "  VCF: ${CURR_VCF}"
    echo "  Truth: ${TRUTH_VCF}"

    # bgzip if needed
    if [[ "${CURR_VCF}" != *.gz ]]; then
        bgzip -c "${CURR_VCF}" > "${CURR_VCF}.gz"
        tabix -f -p vcf "${CURR_VCF}.gz"
        CURR_VCF="${CURR_VCF}.gz"
    fi

    # Build hap.py args
    local HAPPY_ARGS=""

    # For WES mode: restrict to exome targets
    if [[ "${MODE}" == "WES" ]]; then
        HAPPY_ARGS="-T /ref/$(basename ${EXOME_BED})"
    fi

    # Stratification (if available)
    local STRAT_ARGS=""
    if [[ -f "${STRATIFICATION_TSV}" ]]; then
        STRAT_ARGS="--stratification /ref/$(basename ${STRATIFICATION_TSV})"
    fi

    # Run hap.py via Docker
    docker run --rm \
        -v "$(cd "${RESULTS_DIR}" && pwd):/results:ro" \
        -v "$(cd "${REF_DIR}" && pwd):/ref:ro" \
        -v "$(cd "${DATA_DIR}" && pwd):/data:ro" \
        -v "$(cd "${EVAL_OUT}" && pwd):/output" \
        ${HAPPY_IMAGE} \
        /opt/hap.py/bin/hap.py \
        "/data/simulated/${PREFIX}_truth.vcf.gz" \
        "/results/variants/${COV}x/${CALLER_NAME}/$(basename ${CURR_VCF})" \
        -r "/ref/$(basename ${REF_FASTA})" \
        -f "/ref/$(basename ${TRUTH_BED})" \
        ${HAPPY_ARGS} \
        ${STRAT_ARGS} \
        --threads ${THREADS} \
        --engine vcfeval \
        --engine-vcfeval-template "/ref/$(basename ${RTG_SDF})" \
        -o "/output/report"

    echo "  Output: ${EVAL_OUT}/report.*"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}

# ============================================================
# Main: loop over coverages
# ============================================================
mkdir -p "${EVAL_DIR}"

# WGS evaluation — all coverages
for COV in "${COVERAGES_ALL[@]}"; do
    VARIANT_DIR_COV="${RESULTS_DIR}/variants/${COV}x"

    for caller in ${WGS_CALLERS}; do
        VCF=$(get_vcf "${caller}" "${VARIANT_DIR_COV}")
        if [[ -f "${VCF}" ]] || [[ -f "${VCF}.gz" ]]; then
            make_comparison "${VCF}" "${caller}" "${COV}" "WGS"
        else
            echo "SKIP: ${caller} @ ${COV}x WGS - VCF not found: ${VCF}"
        fi
    done
done

# WES evaluation — only high-coverage levels
for COV in "${COVERAGES_WES[@]}" 50; do
    VARIANT_DIR_COV="${RESULTS_DIR}/variants/${COV}x"

    for caller in ${WES_CALLERS}; do
        VCF=$(get_vcf "${caller}" "${VARIANT_DIR_COV}")
        if [[ -f "${VCF}" ]] || [[ -f "${VCF}.gz" ]]; then
            make_comparison "${VCF}" "${caller}" "${COV}" "WES"
        else
            echo "SKIP: ${caller} @ ${COV}x WES - VCF not found: ${VCF}"
        fi
    done
done

echo ""
echo "All evaluations complete."
echo "Results in: ${EVAL_DIR}/"
echo " Structure: ${EVAL_DIR}/{COV}x/{WGS|WES}/{caller}_eval_data/report.*"

# ============================================================
# Gather stats into single TSV for R analysis
# ============================================================
echo ""
echo "Gathering stats..."

ALL_STATS="${EVAL_DIR}/all_stats.tsv"
FIRST=true

for REPORT in $(find "${EVAL_DIR}" -name "report.extended.csv" 2>/dev/null); do
    # Extract coverage and mode from path
    COV=$(echo "${REPORT}" | grep -oP '\d+x' | head -1)
    MODE=$(echo "${REPORT}" | grep -oP '(WGS|WES)' | head -1)
    CALLER=$(basename $(dirname "${REPORT}") | sed 's/_eval_data//')

    if ${FIRST}; then
        # Header + metadata columns
        head -1 "${REPORT}" | sed "s/^/Coverage\tMode\tCaller\t/" > "${ALL_STATS}"
        FIRST=false
    fi

    # Data rows
    tail -n +2 "${REPORT}" | sed "s/^/${COV}\t${MODE}\t${CALLER}\t/" >> "${ALL_STATS}"
done

echo "All stats written to: ${ALL_STATS}"
