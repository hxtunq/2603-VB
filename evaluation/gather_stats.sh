#!/bin/bash
# Gather hap.py evaluation statistics into a single TSV.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"

EVAL_ROOT="${1:-${EVAL_DIR}}"
OUTPUT_TSV="${2:-${EVAL_ROOT}/all_stats.tsv}"

FIRST_REPORT=$(find "${EVAL_ROOT}" -name "report.extended.csv" 2>/dev/null | sort | head -n1)
if [[ -z "${FIRST_REPORT}" ]]; then
    echo "ERROR: No report.extended.csv found in ${EVAL_ROOT}/"
    exit 1
fi

echo "Gathering stats from ${EVAL_ROOT}..."

head -n1 "${FIRST_REPORT}" | sed 's/^/Coverage\tMode\tCaller\t/' | tr ',' '\t' > "${OUTPUT_TSV}"

find "${EVAL_ROOT}" -name "report.extended.csv" | sort | while IFS= read -r report; do
    rel_path="${report#${EVAL_ROOT}/}"
    cov=$(echo "${rel_path}" | cut -d'/' -f1)
    mode=$(echo "${rel_path}" | cut -d'/' -f2)
    caller=$(echo "${rel_path}" | cut -d'/' -f3 | sed 's/_eval_data//')

    tail -n +2 "${report}" | while IFS= read -r line; do
        printf '%s\t%s\t%s\t%s\n' "${cov}" "${mode}" "${caller}" "$(echo "${line}" | tr ',' '\t')"
    done >> "${OUTPUT_TSV}"
done

echo "Stats written to: ${OUTPUT_TSV}"
echo "Total lines: $(wc -l < "${OUTPUT_TSV}")"
