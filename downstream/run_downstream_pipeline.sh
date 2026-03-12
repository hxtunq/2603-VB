#!/usr/bin/env bash
# End-to-end downstream analysis runner.
# This script only invokes downstream scripts and reuses existing outputs as input.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

PYTHON_BIN=${PYTHON_BIN:-python}
R_BIN=${R_BIN:-Rscript}
RUN_ALPHAGENOME=${RUN_ALPHAGENOME:-0}

echo "=========================================="
echo " Downstream Analysis Pipeline"
echo "=========================================="
echo "Project: ${PROJECT_DIR}"
echo "RUN_ALPHAGENOME=${RUN_ALPHAGENOME}"
echo ""

"${PYTHON_BIN}" "${SCRIPT_DIR}/extract_error_patterns.py" --project-root "${PROJECT_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/annotate_variants_strata.py" --project-root "${PROJECT_DIR}"

if [[ "${RUN_ALPHAGENOME}" == "1" ]]; then
  bash "${SCRIPT_DIR}/run_alphagenome_batch.sh"
else
  echo "[INFO] Skip AlphaGenome scoring (set RUN_ALPHAGENOME=1 to enable)"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/aggregate_alphagenome_scores.py" --project-root "${PROJECT_DIR}"

"${R_BIN}" "${SCRIPT_DIR}/upset_plot_enhanced.R" "${PROJECT_DIR}"
"${R_BIN}" "${SCRIPT_DIR}/plot_strata_characterization.R" "${PROJECT_DIR}"
"${R_BIN}" "${SCRIPT_DIR}/analyze_coverage_sensitivity.R" "${PROJECT_DIR}"

echo ""
echo "[DONE] Downstream pipeline finished"
