#!/bin/bash
#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $1" >&2
}

ensure_dir() {
    mkdir -p "$1"
}

start_timer() {
    export STEP_START_TIME=$(date +%s)
}

end_timer() {
    local step_name="$1"
    local end_time=$(date +%s)
    local duration=$((end_time - STEP_START_TIME))
    log_info "${step_name} completed in ${duration} seconds"
    echo "${step_name},${duration}" >> "${LOG_DIR}/runtime.csv"
}

parse_time_to_seconds() {
    local wall_time="$1"
    local h=0
    local m=0
    local s=0

    if [[ "${wall_time}" == *:*:* ]]; then
        IFS=':' read -r h m s <<< "${wall_time}"
    elif [[ "${wall_time}" == *:* ]]; then
        IFS=':' read -r m s <<< "${wall_time}"
    else
        s="${wall_time}"
    fi

    awk -v h="${h}" -v m="${m}" -v s="${s}" 'BEGIN {printf "%.3f", (h*3600)+(m*60)+s}'
}

append_resource_metrics() {
    local caller="$1"
    local step_name="$2"
    local time_log="$3"
    local metrics_file="${LOG_DIR}/resource_usage.tsv"
    local header="Caller\tStep\tWallClockSeconds\tCPUPercent\tMaxRSS_GB"

    local wall_raw
    local cpu_raw
    local rss_kb

    wall_raw=$(awk -F': *' '/Elapsed \(wall clock\) time/{print $2; exit}' "${time_log}" || true)
    cpu_raw=$(awk -F': *' '/Percent of CPU this job got/{print $2; exit}' "${time_log}" | tr -d '%' || true)
    rss_kb=$(awk -F': *' '/Maximum resident set size \(kbytes\)/{print $2; exit}' "${time_log}" || true)

    local wall_seconds="NA"
    local cpu_percent="NA"
    local rss_gb="NA"

    if [[ -n "${wall_raw}" ]]; then
        wall_seconds=$(parse_time_to_seconds "${wall_raw}")
    fi
    if [[ -n "${cpu_raw}" ]]; then
        cpu_percent="${cpu_raw}"
    fi
    if [[ -n "${rss_kb}" ]]; then
        rss_gb=$(awk -v kb="${rss_kb}" 'BEGIN {printf "%.3f", kb/1024/1024}')
        local rss_limit_kb=$((14 * 1024 * 1024))
        if [[ "${rss_kb}" -gt "${rss_limit_kb}" ]]; then
            log_warn "${caller}/${step_name} exceeded RAM cap 14G (MaxRSS=${rss_gb}G)"
        fi
    fi

    ensure_dir "$(dirname "${metrics_file}")"
    {
        echo -e "${header}"
        awk 'NR>1 {print}' "${metrics_file}" 2>/dev/null || true
        echo -e "${caller}\t${step_name}\t${wall_seconds}\t${cpu_percent}\t${rss_gb}"
    } > "${metrics_file}.tmp"
    mv "${metrics_file}.tmp" "${metrics_file}"
}
