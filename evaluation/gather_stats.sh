#!/bin/bash
# Gather hap.py evaluation statistics into a single TSV
# Based on: caller_benchmark-main/gather_stats.sh

set -e

source "$(dirname "$0")/../config/config.sh"

OUTPUT_TSV="${EVAL_DIR}/all_stats.tsv"

# Print header from first report
FIRST_REPORT=$(find "${EVAL_DIR}" -name "report.extended.csv" | head -n1)
if [[ -z "${FIRST_REPORT}" ]]; then
    echo "ERROR: No report.extended.csv found in ${EVAL_DIR}/"
    exit 1
fi

echo "Gathering stats from ${EVAL_DIR}..."

# Header: prepend Caller column
head -n1 "${FIRST_REPORT}" | sed 's/^/Caller\t/' | tr ',' '\t' > "${OUTPUT_TSV}"

# Data rows
for eval_dir in "${EVAL_DIR}"/*_eval_data; do
    if [[ ! -d "${eval_dir}" ]]; then continue; fi

    CALLER=$(basename "${eval_dir}" | sed 's/_eval_data//')
    REPORT="${eval_dir}/report.extended.csv"

    if [[ -f "${REPORT}" ]]; then
        tail -n +2 "${REPORT}" | while IFS= read -r line; do
            echo -e "${CALLER}\t$(echo "${line}" | tr ',' '\t')"
        done >> "${OUTPUT_TSV}"
    else
        echo "WARN: No report found for ${CALLER}"
    fi
done

echo "Stats written to: ${OUTPUT_TSV}"
echo "Total lines: $(wc -l < "${OUTPUT_TSV}")"
