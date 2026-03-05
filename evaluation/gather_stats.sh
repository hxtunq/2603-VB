#!/bin/bash
# Gather hap.py evaluation statistics into a single TSV
# Based on: caller_benchmark-main/gather_stats.sh
# Supports multi-coverage structure: ${EVAL_DIR}/{COV}x/{WGS|WES}/{caller}_eval_data/

set -e

source "$(dirname "$0")/../config/config.sh"

OUTPUT_TSV="${EVAL_DIR}/all_stats.tsv"

# Find first report for header
FIRST_REPORT=$(find "${EVAL_DIR}" -name "report.extended.csv" 2>/dev/null | head -n1)
if [[ -z "${FIRST_REPORT}" ]]; then
    echo "ERROR: No report.extended.csv found in ${EVAL_DIR}/"
    exit 1
fi

echo "Gathering stats from ${EVAL_DIR}..."

# Header: prepend Coverage, Mode, Caller columns
head -n1 "${FIRST_REPORT}" | sed 's/^/Coverage\tMode\tCaller\t/' | tr ',' '\t' > "${OUTPUT_TSV}"

# Data rows — walk all report files recursively
find "${EVAL_DIR}" -name "report.extended.csv" | sort | while IFS= read -r REPORT; do
    # Extract coverage, mode, caller from path
    # Expected path: ${EVAL_DIR}/{COV}x({_wes})/{WGS|WES}/{caller}_eval_data/report.extended.csv
    REL_PATH="${REPORT#${EVAL_DIR}/}"
    COV=$(echo "${REL_PATH}" | cut -d'/' -f1)         # e.g. 10x or 50x_wes
    MODE=$(echo "${REL_PATH}" | cut -d'/' -f2)         # e.g. WGS or WES
    CALLER=$(echo "${REL_PATH}" | cut -d'/' -f3 | sed 's/_eval_data//')  # e.g. gatk

    tail -n +2 "${REPORT}" | while IFS= read -r line; do
        echo -e "${COV}\t${MODE}\t${CALLER}\t$(echo "${line}" | tr ',' '\t')"
    done >> "${OUTPUT_TSV}"
done

echo "Stats written to: ${OUTPUT_TSV}"
echo "Total lines: $(wc -l < "${OUTPUT_TSV}")"
