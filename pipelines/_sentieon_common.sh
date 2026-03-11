#!/bin/bash

sentieon_die() {
    echo "ERROR: $*" >&2
    exit 1
}

sentieon_skip_if_no_license() {
    local label="$1"
    if [[ -z "${SENTIEON_LICENSE:-}" ]]; then
        echo "WARNING: SENTIEON_LICENSE not set. Skipping ${label}."
        exit 0
    fi
}

sentieon_require_command() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || sentieon_die "'${cmd}' not found in PATH"
}

sentieon_require_file() {
    local path="$1"
    [[ -f "${path}" ]] || sentieon_die "Required file not found: ${path}"
}

sentieon_require_vcf_index() {
    local path="$1"

    case "${path}" in
        *.vcf.gz)
            [[ -f "${path}.tbi" || -f "${path}.csi" ]] || \
                sentieon_die "Required VCF index not found for: ${path}"
            ;;
        *.vcf)
            [[ -f "${path}.idx" || -f "${path}.tbi" || -f "${path}.csi" ]] || \
                sentieon_die "Required VCF index not found for: ${path}"
            ;;
    esac
}

sentieon_require_reference_fasta() {
    local fasta="$1"

    sentieon_require_file "${fasta}"
    sentieon_require_file "${fasta}.fai"
}

sentieon_require_bwa_index() {
    local fasta="$1"

    for ext in amb ann bwt pac sa; do
        sentieon_require_file "${fasta}.${ext}"
    done
}

sentieon_resolve_model_bundle() {
    local input="$1"
    local need_bwa="${2:-false}"
    local parent=""
    local base=""
    local candidate=""

    if [[ -f "${input}" ]]; then
        printf '%s\n' "${input}"
        return
    fi

    if [[ -d "${input}" ]]; then
        parent=$(cd "$(dirname "${input}")" >/dev/null 2>&1 && pwd)
        base=$(basename "${input}")

        for candidate in \
            "${input}.bundle" \
            "${parent}/${base}.bundle" \
            "${parent}/${base%.bundle}.bundle"
        do
            if [[ -f "${candidate}" ]]; then
                echo "WARNING: '${input}' is an unpacked Sentieon model directory. Using sibling bundle archive '${candidate}' for sentieon-cli." >&2
                printf '%s\n' "${candidate}"
                return
            fi
        done

        if [[ -f "${input}/dnascope.model" ]]; then
            if [[ "${need_bwa}" == "true" ]]; then
                sentieon_require_file "${input}/bwa.model"
            fi
            sentieon_die "DNAscope model path must point to a .bundle file for sentieon-cli, not an unpacked directory: ${input}"
        fi

        sentieon_die "DNAscope model bundle path is a directory, but sentieon-cli expects a .bundle file: ${input}"
    fi

    sentieon_die "Required Sentieon model bundle file not found: ${input}"
}

sentieon_require_dnascope_cli_stack() {
    sentieon_require_command sentieon-cli
    if ! command -v sentieon >/dev/null 2>&1; then
        sentieon_die "'sentieon' not found in PATH. sentieon-cli dnascope shells out to 'sentieon driver'. Add the Sentieon genomics bin directory to PATH or set SENTIEON_BIN_DIR before sourcing config/config.sh."
    fi
}

sentieon_prepare_layout() {
    local cov_suffix="$1"
    local caller="$2"

    OUT_DIR="${VARIANT_DIR}/${cov_suffix}/${caller}"
    TIMEDIR="${LOG_DIR}/${cov_suffix}/time"
    METRICS="${LOG_DIR}/${cov_suffix}/benchmark_metrics.tsv"

    mkdir -p "${OUT_DIR}" "${TIMEDIR}"

    if [[ ! -f "${METRICS}" ]]; then
        printf 'Caller\tPipeline\tWallClock_sec\tCPU_percent\tMaxRSS_kB\n' > "${METRICS}"
    fi
}

sentieon_log_metrics() {
    local caller="$1"
    local pipeline="$2"
    local timefile="$3"
    local wall cpu rss secs

    sentieon_require_file "${timefile}"

    wall=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "${timefile}")
    cpu=$(awk -F': ' '/Percent of CPU this job got/{print $2}' "${timefile}" | tr -d '%')
    rss=$(awk -F': ' '/Maximum resident set size \(kbytes\)/{print $2}' "${timefile}")
    secs=$(echo "${wall}" | awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}')

    printf '%s\t%s\t%s\t%s\t%s\n' "${caller}" "${pipeline}" "${secs}" "${cpu}" "${rss}" >> "${METRICS}"
    echo "  [METRICS] ${pipeline}: ${wall} wall, ${cpu}% CPU, MaxRSS=${rss} kB"
}

sentieon_set_fastq_paths() {
    local cov="$1"

    FASTQ_R1="${SIM_DIR}/${PREFIX}_${cov}x_R1.fastq.gz"
    FASTQ_R2="${SIM_DIR}/${PREFIX}_${cov}x_R2.fastq.gz"
    COV_SUFFIX="${cov}x"
}

sentieon_build_readgroup() {
    local cov="$1"
    local rg_id="${SAMPLE_NAME}_${cov}x"

    printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
        "${rg_id}" "${SAMPLE_NAME}"
}