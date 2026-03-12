#!/usr/bin/env bash
# Run AlphaGenome scoring for selected callers, coverages, and FN/FP types.
# This wrapper does not modify existing scripts and only orchestrates calls
# to batch_alphagenome_fn_fp.py.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

PYTHON_BIN=${PYTHON_BIN:-python}
COVERAGES=${COVERAGES:-"10 20 30 50"}
CALLERS=${CALLERS:-"gatk deepvariant strelka2 freebayes dnascope"}
ERROR_TYPES=${ERROR_TYPES:-"FN FP"}
SEQUENCE_LENGTH=${SEQUENCE_LENGTH:-"1MB"}
SCORERS=${SCORERS:-"all"}
BATCH_SIZE=${BATCH_SIZE:-50}
SKIP_IF_EXISTS=${SKIP_IF_EXISTS:-1}
PATTERN_PREFIXES=${PATTERN_PREFIXES:-"all_5,single_caller"}
PATTERNS_DIR=${PATTERNS_DIR:-"results/analysis/error_patterns"}

OUT_DIR="${PROJECT_DIR}/results/alphagenome"
mkdir -p "${OUT_DIR}"

if [[ -z "${ALPHAGENOME_API_KEY:-}" ]]; then
  echo "[WARN] ALPHAGENOME_API_KEY is not set. The downstream scorer may fail."
fi

echo "=========================================="
echo " Downstream AlphaGenome Batch Runner"
echo "=========================================="
echo "Project: ${PROJECT_DIR}"
echo "Callers: ${CALLERS}"
echo "Coverages: ${COVERAGES}"
echo "Error types: ${ERROR_TYPES}"
echo "Pattern prefixes: ${PATTERN_PREFIXES:-<none>}"
echo "Patterns dir: ${PATTERNS_DIR}"
echo "Seq length: ${SEQUENCE_LENGTH}"
echo "Batch size: ${BATCH_SIZE}"
echo ""

for caller in ${CALLERS}; do
  for cov in ${COVERAGES}; do
    for err in ${ERROR_TYPES}; do
      out_csv="${OUT_DIR}/${caller}_${cov}x_${err}_variant_scores_tidy.csv"

      if [[ "${SKIP_IF_EXISTS}" == "1" && -f "${out_csv}" ]]; then
        echo "[SKIP] ${out_csv} already exists"
        continue
      fi

      echo "[RUN] caller=${caller} cov=${cov}x err=${err}"
      "${PYTHON_BIN}" "${PROJECT_DIR}/batch_alphagenome_fn_fp.py" \
        --caller "${caller}" \
        --error-type "${err}" \
        --coverage "${cov}" \
        --patterns-dir "${PATTERNS_DIR}" \
        --pattern-prefixes "${PATTERN_PREFIXES}" \
        --sequence-length "${SEQUENCE_LENGTH}" \
        --scorers "${SCORERS}" \
        --batch-size "${BATCH_SIZE}" \
        --root "${PROJECT_DIR}"
    done
  done
done

echo ""
echo "[DONE] AlphaGenome batch orchestration finished"
