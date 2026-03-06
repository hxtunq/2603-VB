#!/usr/bin/env bash

happy_die() {
    echo "ERROR: $*" >&2
    exit 1
}

happy_require_command() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || happy_die "'${cmd}' not found in PATH"
}

happy_require_file() {
    local path="$1"
    [[ -e "${path}" ]] || happy_die "Required path not found: ${path}"
}

happy_prepare_query_vcf() {
    local curr_vcf="$1"

    if [[ "${curr_vcf}" != *.gz ]]; then
        bgzip -c "${curr_vcf}" > "${curr_vcf}.gz"
        tabix -f -p vcf "${curr_vcf}.gz"
        curr_vcf="${curr_vcf}.gz"
    fi

    printf '%s\n' "${curr_vcf}"
}

happy_cov_suffix() {
    local cov="$1"
    local mode="$2"

    if [[ "${mode}" == "WES" ]]; then
        printf '%sx_wes\n' "${cov}"
    else
        printf '%sx\n' "${cov}"
    fi
}

happy_make_comparison() {
    local curr_vcf="$1"
    local caller_name="$2"
    local cov="$3"
    local mode="$4"
    local eval_root="$5"
    local cov_suffix eval_out
    local happy_args=()
    local strat_args=()

    cov_suffix=$(happy_cov_suffix "${cov}" "${mode}")
    eval_out="${eval_root}/${cov_suffix}/${mode}/${caller_name}_eval_data"
    mkdir -p "${eval_out}"

    curr_vcf=$(happy_prepare_query_vcf "${curr_vcf}")

    if [[ "${mode}" == "WES" ]]; then
        happy_args=(-T "/ref/$(basename "${EXOME_BED}")")
    fi

    if [[ -f "${STRATIFICATION_TSV}" ]]; then
        strat_args=(--stratification "/ref/$(basename "${STRATIFICATION_TSV}")")
    fi

    echo "Evaluating: ${caller_name} @ ${cov}x (${mode})"
    echo "  VCF: ${curr_vcf}"
    echo "  Truth: ${TRUTH_VCF}"

    docker run --rm \
        -v "$(cd "${RESULTS_DIR}" && pwd):/results:ro" \
        -v "$(cd "${REF_DIR}" && pwd):/ref:ro" \
        -v "$(cd "${DATA_DIR}" && pwd):/data:ro" \
        -v "$(cd "${eval_out}" && pwd):/output" \
        "${HAPPY_IMAGE}" \
        /opt/hap.py/bin/hap.py \
        "/data/simulated/${PREFIX}_truth.vcf.gz" \
        "/results/variants/${cov_suffix}/${caller_name}/$(basename "${curr_vcf}")" \
        -r "/ref/$(basename "${REF_FASTA}")" \
        -f "/ref/$(basename "${TRUTH_BED}")" \
        "${happy_args[@]}" \
        "${strat_args[@]}" \
        --threads "${THREADS}" \
        --engine vcfeval \
        --engine-vcfeval-template "/ref/$(basename "${RTG_SDF}")" \
        --preprocess-truth \
        --no-fixchr \
        -o "/output/report"

    echo "  Output: ${eval_out}/report.*"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}
