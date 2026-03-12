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

    # --- PASS filter + Normalize for fair evaluation ---
    local vcf_dir vcf_base norm_vcf
    vcf_dir=$(dirname "${curr_vcf}")
    vcf_base=$(basename "${curr_vcf}" .vcf.gz)
    norm_vcf="${vcf_dir}/${vcf_base}.norm.vcf.gz"

    if [[ ! -f "${norm_vcf}" ]]; then
        echo "  [PREP] Filtering PASS + normalizing: $(basename "${curr_vcf}")"
        bcftools view -f 'PASS,.' "${curr_vcf}" \
            | bcftools norm -m -both -f "${REF_FASTA}" \
            | bgzip -c > "${norm_vcf}"
        tabix -f -p vcf "${norm_vcf}"
    fi
    curr_vcf="${norm_vcf}"

    printf '%s\n' "${curr_vcf}"
}

happy_prepare_truth_vcf() {
    local truth_prepared="${EVAL_DIR}/truth_prepared.vcf.gz"
    local truth_tmp="${truth_prepared}.tmp"
    local sample_count=0

    # Ensure EVAL_DIR exists before making subdirectories
    mkdir -p "$(dirname "${truth_prepared}")"

    if [[ -s "${truth_prepared}" && -s "${truth_prepared}.tbi" && "${truth_prepared}" -nt "${TRUTH_VCF}" ]]; then
        printf '%s\n' "${truth_prepared}"
        return
    fi

    rm -f "${truth_tmp}" "${truth_tmp}.tbi"

    sample_count=$(bcftools query -l "${TRUTH_VCF}" | awk 'NF {c++} END {print c+0}')

    if (( sample_count > 0 )); then
        echo "  [TRUTH] Input truth VCF has ${sample_count} sample(s); keeping only the first sample column." >&2
        zcat "${TRUTH_VCF}" | awk 'BEGIN{OFS="\t"; ff=0; gt=0}
          /^##fileformat=/ {ff=1; print; next}
          /^##FORMAT=<ID=GT,/ {gt=1; print; next}
          /^##/ {print; next}
          /^#CHROM/ {
            if(ff==0) print "##fileformat=VCFv4.2"
            if(gt==0) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            print $1,$2,$3,$4,$5,$6,$7,$8,"FORMAT","TRUTH"
            next
          }
          {
            if(NF < 10) {
              printf "ERROR: expected FORMAT + sample columns at %s:%s but found %d fields\n", $1, $2, NF > "/dev/stderr"
              exit 1
            }
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10
          }
        ' | bcftools norm -f "${REF_FASTA}" -m -both \
          | bgzip -c > "${truth_tmp}"
    else
        echo "  [TRUTH] Input truth VCF is sites-only; adding GT=1/1 sample column." >&2
        zcat "${TRUTH_VCF}" | awk 'BEGIN{OFS="\t"; ff=0}
          /^##fileformat=/ {ff=1; print; next}
          /^##/ {print; next}
          /^#CHROM/ {
            if(ff==0) print "##fileformat=VCFv4.2"
            print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            print $0,"FORMAT","TRUTH"
            next
          }
          {print $0,"GT","1/1"}
        ' | bcftools norm -f "${REF_FASTA}" -m -both \
          | bgzip -c > "${truth_tmp}"
    fi

    tabix -f -p vcf "${truth_tmp}"
    mv -f "${truth_tmp}" "${truth_prepared}"
    mv -f "${truth_tmp}.tbi" "${truth_prepared}.tbi"
    
    printf '%s\n' "${truth_prepared}"
}

happy_cov_suffix() {
    local cov="$1"
    printf '%sx\n' "${cov}"
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

    cov_suffix=$(happy_cov_suffix "${cov}")
    eval_out="${eval_root}/${cov_suffix}/${mode}/${caller_name}_eval_data"
    mkdir -p "${eval_out}"

    curr_vcf=$(happy_prepare_query_vcf "${curr_vcf}")

    local truth_vcf
    truth_vcf=$(happy_prepare_truth_vcf)

    if [[ -f "${STRATIFICATION_TSV}" ]]; then
        strat_args=(--stratification "/ref/$(basename "${STRATIFICATION_TSV}")")
    fi

    echo "Evaluating: ${caller_name} @ ${cov}x (${mode})"
    echo "  VCF: ${curr_vcf}"
    echo "  Truth: ${truth_vcf}"

    docker run --rm \
        -v "$(cd "${RESULTS_DIR}" && pwd):/results:ro" \
        -v "$(cd "${REF_DIR}" && pwd):/ref:ro" \
        -v "$(cd "${DATA_DIR}" && pwd):/data:ro" \
        -v "$(cd "${eval_out}" && pwd):/output" \
        "${HAPPY_IMAGE}" \
        /opt/hap.py/bin/hap.py \
        "/results/eval/$(basename "${truth_vcf}")" \
        "/results/variants/${cov_suffix}/${caller_name}/$(basename "${curr_vcf}")" \
        -r "/ref/$(basename "${REF_FASTA}")" \
        -f "/ref/$(basename "${TRUTH_BED}")" \
        "${happy_args[@]}" \
        "${strat_args[@]}" \
        --threads "${THREADS}" \
        --engine vcfeval \
        --engine-vcfeval-template "/ref/$(basename "${RTG_SDF}")" \
        -o "/output/report"

    echo "  Output: ${eval_out}/report.*"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
}
