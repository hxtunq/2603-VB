#!/usr/bin/env bash
# End-to-end downstream analysis runner.
# This script only invokes downstream scripts and reuses existing outputs as input.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

PYTHON_BIN=${PYTHON_BIN:-python}
R_BIN=${R_BIN:-Rscript}
RUN_ALPHAGENOME=${RUN_ALPHAGENOME:-0}
CONCORDANCE_DIR=${CONCORDANCE_DIR:-results/eval/vcfeval}
CONCORDANCE_SEARCH_ROOT=${CONCORDANCE_SEARCH_ROOT:-results/eval}
AUTO_BUILD_CONCORDANCE=${AUTO_BUILD_CONCORDANCE:-1}

echo "=========================================="
echo " Downstream Analysis Pipeline"
echo "=========================================="
echo "Project: ${PROJECT_DIR}"
echo "RUN_ALPHAGENOME=${RUN_ALPHAGENOME}"
echo "CONCORDANCE_DIR=${CONCORDANCE_DIR}"
echo "CONCORDANCE_SEARCH_ROOT=${CONCORDANCE_SEARCH_ROOT}"
echo "AUTO_BUILD_CONCORDANCE=${AUTO_BUILD_CONCORDANCE}"
echo ""

if [[ "${AUTO_BUILD_CONCORDANCE}" == "1" ]]; then
  "${PYTHON_BIN}" "${SCRIPT_DIR}/build_concordance_from_vcfs.py" \
    --project-root "${PROJECT_DIR}" \
    --search-root "${CONCORDANCE_SEARCH_ROOT}" \
    --out-dir "${CONCORDANCE_DIR}"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/extract_error_patterns.py" --project-root "${PROJECT_DIR}" --vcfeval-dir "${CONCORDANCE_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/annotate_variants_strata.py" --project-root "${PROJECT_DIR}"

if [[ "${RUN_ALPHAGENOME}" == "1" ]]; then
  bash "${SCRIPT_DIR}/run_alphagenome_batch.sh"
else
  echo "[INFO] Skip AlphaGenome scoring (set RUN_ALPHAGENOME=1 to enable)"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/aggregate_alphagenome_scores.py" --project-root "${PROJECT_DIR}"

"${R_BIN}" "${SCRIPT_DIR}/upset_plot_enhanced.R" "${PROJECT_DIR}" "${PROJECT_DIR}/${CONCORDANCE_DIR}"
"${R_BIN}" "${SCRIPT_DIR}/plot_strata_characterization.R" "${PROJECT_DIR}"
"${R_BIN}" "${SCRIPT_DIR}/analyze_coverage_sensitivity.R" "${PROJECT_DIR}" "${PROJECT_DIR}/${CONCORDANCE_DIR}"

echo ""
echo "[DONE] Downstream pipeline finished"
