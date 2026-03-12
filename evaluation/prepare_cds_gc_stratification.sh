#!/usr/bin/env bash
# =============================================================================
# Prepare CDS-specific GC stratification BED files.
#
# Generates:
#   - CDS_GC_{BIN}.bed by computing %GC inside the CDS canonical regions.
#   - Appends these BED paths to stratification_chr22.tsv.
#
# Usage:
#   bash evaluation/prepare_cds_gc_stratification.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

STRAT_DIR="${REF_DIR}/stratification"
GC_DIR="${REF_DIR}/gc_strata"
mkdir -p "${GC_DIR}"

CDS_CANONICAL_BED="${STRAT_DIR}/CDS-canonical.${CHR_TO_USE}.bed"

if [[ ! -f "${CDS_CANONICAL_BED}" ]]; then
    echo "ERROR: CDS canonical BED not found at ${CDS_CANONICAL_BED}"
    echo "Please run evaluation/prepare_stratification_hg38.sh first."
    exit 1
fi

echo "=========================================="
echo " Preparing %GC in CDS Region Strata"
echo "=========================================="

echo "Running generate_cds_gc_strata.py..."
python3 "${SCRIPT_DIR}/generate_cds_gc_strata.py" "${REF_FASTA}" "${CDS_CANONICAL_BED}" "${GC_DIR}"

echo ""
echo "=========================================="
echo " CDS GC Stratification complete"
echo "=========================================="
echo "Next: Re-run bash evaluation/eval_happy.sh to apply the new strata"
