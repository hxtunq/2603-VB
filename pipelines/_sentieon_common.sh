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

sentieon_require_dir() {
    local path="$1"
    [[ -d "${path}" ]] || sentieon_die "Required directory not found: ${path}"
}

sentieon_require_prefix() {
    local prefix="$1"
    compgen -G "${prefix}*" >/dev/null || sentieon_die "No files found for index prefix: ${prefix}"
}

sentieon_abs_dir() {
    local path="$1"
    (cd "${path}" >/dev/null 2>&1 && pwd)
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

sentieon_init_license_args() {
    SENTIEON_DOCKER_LICENSE_ARGS=()

    if [[ -f "${SENTIEON_LICENSE}" ]]; then
        local license_dir license_name
        license_dir=$(sentieon_abs_dir "$(dirname "${SENTIEON_LICENSE}")")
        license_name=$(basename "${SENTIEON_LICENSE}")
        SENTIEON_DOCKER_LICENSE_ARGS=(
            -e "SENTIEON_LICENSE=/licenses/${license_name}"
            -v "${license_dir}:/licenses:ro"
        )
    else
        SENTIEON_DOCKER_LICENSE_ARGS=(-e "SENTIEON_LICENSE=${SENTIEON_LICENSE}")
    fi
}

sentieon_set_pcr_arg() {
    SENTIEON_PCR_ARG=""
    if [[ "${PCRFREE:-false}" == "true" ]]; then
        SENTIEON_PCR_ARG="--pcr_indel_model none"
    fi
}

sentieon_set_fastq_paths() {
    local cov="$1"
    local mode="$2"

    if [[ "${mode}" == "wes" ]]; then
        FASTQ_R1="${SIM_DIR}/${PREFIX}_${cov}x_wes_R1.fastq.gz"
        FASTQ_R2="${SIM_DIR}/${PREFIX}_${cov}x_wes_R2.fastq.gz"
        COV_SUFFIX="${cov}x_wes"
    else
        FASTQ_R1="${SIM_DIR}/${PREFIX}_${cov}x_R1.fastq.gz"
        FASTQ_R2="${SIM_DIR}/${PREFIX}_${cov}x_R2.fastq.gz"
        COV_SUFFIX="${cov}x"
    fi
}
