#!/usr/bin/env bash
# =============================================================================
# Prepare stratification BED files for hg38 / chr22 benchmarking.
#
# Downloads GENCODE v44 GTF, generates:
#   - CDS (protein-coding exon) BED files with 0/25/50/100bp padding
#   - Coverage-stratification BEDs from BAM depth
#   - Assembles stratification_chr22.tsv for hap.py
#
# GC strata are handled separately by generate_gc_strata.py.
#
# Usage:
#   bash evaluation/prepare_stratification_hg38.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

STRAT_DIR="${REF_DIR}/stratification"
mkdir -p "${STRAT_DIR}"

GENCODE_VERSION="44"
GENCODE_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
GENCODE_GTF="${STRAT_DIR}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"

echo "=========================================="
echo " Preparing Stratification BEDs (hg38)"
echo "=========================================="

# =============================================================================
# 1) GENCODE CDS regions
# =============================================================================

echo ""
echo "--- Step 1: GENCODE CDS regions ---"

if [[ ! -f "${GENCODE_GTF}" ]]; then
    echo "  Downloading GENCODE v${GENCODE_VERSION} GTF..."
    wget -q -O "${GENCODE_GTF}" "${GENCODE_GTF_URL}"
    echo "  Downloaded: ${GENCODE_GTF}"
else
    echo "  Already exists: ${GENCODE_GTF}"
fi

# Extract protein-coding CDS on chr22 only
CDS_BED="${STRAT_DIR}/${CHR_TO_USE}_proteincoding_only.bed"
if [[ ! -f "${CDS_BED}" ]]; then
    echo "  Extracting protein-coding CDS for ${CHR_TO_USE}..."
    zcat "${GENCODE_GTF}" \
        | awk -v chr="${CHR_TO_USE}" '$1 == chr && $3 == "CDS" && $0 ~ /gene_type "protein_coding"/' \
        | awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "${CDS_BED}"
    echo "  Created: ${CDS_BED} ($(wc -l < "${CDS_BED}") intervals)"
else
    echo "  Already exists: ${CDS_BED}"
fi

# Extract canonical protein-coding CDS on chr22 only
CDS_CANONICAL_BED="${STRAT_DIR}/CDS-canonical.${CHR_TO_USE}.bed"
if [[ ! -f "${CDS_CANONICAL_BED}" ]]; then
    echo "  Extracting canonical protein-coding CDS for ${CHR_TO_USE}..."
    zcat "${GENCODE_GTF}" \
        | awk -v chr="${CHR_TO_USE}" '$1 == chr && $3 == "CDS" && $0 ~ /gene_type "protein_coding"/ && $0 ~ /tag "Ensembl_canonical"/' \
        | awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "${CDS_CANONICAL_BED}"
    echo "  Created: ${CDS_CANONICAL_BED} ($(wc -l < "${CDS_CANONICAL_BED}") intervals)"
else
    echo "  Already exists: ${CDS_CANONICAL_BED}"
fi

# Generate CDS vicinity (padding) BED files
for PAD in 0 25 50 100; do
    VICINITY_BED="${STRAT_DIR}/${CHR_TO_USE}_cds_vicinity_${PAD}bp.bed"
    if [[ ! -f "${VICINITY_BED}" ]]; then
        echo "  Generating CDS ± ${PAD}bp..."
        if [[ "${PAD}" -eq 0 ]]; then
            cp "${CDS_BED}" "${VICINITY_BED}"
        else
            bedtools slop -i "${CDS_BED}" -g "${REF_FAI}" -b "${PAD}" \
                | sort -k1,1 -k2,2n \
                | bedtools merge -i - \
                > "${VICINITY_BED}"
        fi
        echo "  Created: ${VICINITY_BED} ($(wc -l < "${VICINITY_BED}") intervals)"
    else
        echo "  Already exists: ${VICINITY_BED}"
    fi
done

# =============================================================================
# 2) Coverage stratification BEDs
# =============================================================================

echo ""
echo "--- Step 2: Coverage stratification BEDs ---"

# Coverage bins (normalized: depth / mean_depth)
COV_BINS=("0_0.25" "0.25_0.5" "0.5_0.75" "0.75_1" "1_1.25" "1.25_1.5" "1.5_1.75" "1.75_2" "2_Inf")

for COV in "${COVERAGES_WGS[@]}"; do
    BAM="${PREPROC_DIR}/${COV}x/${SAMPLE_NAME}_${CHR_TO_USE}_${COV}x_recal.bam"

    if [[ ! -f "${BAM}" ]]; then
        echo "  SKIP: BAM not found for ${COV}x: ${BAM}"
        continue
    fi

    COV_STRAT_DIR="${STRAT_DIR}/coverage_${COV}x"
    if [[ -d "${COV_STRAT_DIR}" ]] && [[ $(ls "${COV_STRAT_DIR}"/*.bed 2>/dev/null | wc -l) -ge 5 ]]; then
        echo "  Coverage strata already exist for ${COV}x — skipping"
        continue
    fi

    mkdir -p "${COV_STRAT_DIR}"

    echo "  Computing depth for ${COV}x BAM..."
    DEPTH_FILE="${COV_STRAT_DIR}/depth.tsv"
    samtools depth -a -r "${CHR_TO_USE}" "${BAM}" > "${DEPTH_FILE}"

    # Compute mean depth
    MEAN_DEPTH=$(awk '{s+=$3; n++} END {printf "%.2f", s/n}' "${DEPTH_FILE}")
    echo "  Mean depth: ${MEAN_DEPTH}"

    # Window-based coverage stratification (1000bp windows)
    WINDOW_SIZE=1000
    echo "  Generating coverage strata (window=${WINDOW_SIZE}bp)..."

    python3 - "${DEPTH_FILE}" "${COV_STRAT_DIR}" "${MEAN_DEPTH}" "${WINDOW_SIZE}" "${CHR_TO_USE}" <<'PYEOF'
import sys, os

depth_file = sys.argv[1]
out_dir    = sys.argv[2]
mean_depth = float(sys.argv[3])
window     = int(sys.argv[4])
chrom      = sys.argv[5]

bins = [
    ("0_0.25",    0.0,  0.25),
    ("0.25_0.5",  0.25, 0.5),
    ("0.5_0.75",  0.5,  0.75),
    ("0.75_1",    0.75, 1.0),
    ("1_1.25",    1.0,  1.25),
    ("1.25_1.5",  1.25, 1.5),
    ("1.5_1.75",  1.5,  1.75),
    ("1.75_2",    1.75, 2.0),
    ("2_Inf",     2.0,  999.0),
]

# Read depth and compute windowed averages
windows = {}
with open(depth_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        pos = int(parts[1])
        depth = int(parts[2])
        w_start = (pos // window) * window
        if w_start not in windows:
            windows[w_start] = [0, 0]
        windows[w_start][0] += depth
        windows[w_start][1] += 1

# Write BED files per bin
handles = {}
for name, _, _ in bins:
    handles[name] = open(os.path.join(out_dir, f"wgs_mean_cov_{name}.bed"), 'w')

for w_start in sorted(windows.keys()):
    total_d, count = windows[w_start]
    avg_depth = total_d / count
    norm = avg_depth / mean_depth if mean_depth > 0 else 0
    w_end = w_start + window

    for name, lo, hi in bins:
        if lo <= norm < hi:
            handles[name].write(f"{chrom}\t{w_start}\t{w_end}\n")
            break

for fh in handles.values():
    fh.close()

for name, _, _ in bins:
    path = os.path.join(out_dir, f"wgs_mean_cov_{name}.bed")
    n = sum(1 for _ in open(path))
    print(f"    {name}: {n} windows")
PYEOF

    rm -f "${DEPTH_FILE}"
done

# =============================================================================
# 3) Assemble stratification TSV
# =============================================================================

echo ""
echo "--- Step 3: Assembling stratification TSV ---"

STRAT_TSV="${REF_DIR}/stratification_chr22.tsv"

# Start fresh
> "${STRAT_TSV}"

# CDS and vicinity
echo "proteincoding_only	stratification/${CHR_TO_USE}_proteincoding_only.bed" >> "${STRAT_TSV}"
echo "CDS_canonical	stratification/CDS-canonical.${CHR_TO_USE}.bed" >> "${STRAT_TSV}"
for PAD in 0 25 50 100; do
    echo "cds_vicinity_${PAD}bp	stratification/${CHR_TO_USE}_cds_vicinity_${PAD}bp.bed" >> "${STRAT_TSV}"
done

# GC strata (if they exist from generate_gc_strata.py)
GC_DIR="${REF_DIR}/gc_strata"
if [[ -d "${GC_DIR}" ]]; then
    for gc_bed in "${GC_DIR}"/*.bed; do
        name=$(basename "${gc_bed}" .bed)
        echo "${name}	gc_strata/$(basename "${gc_bed}")" >> "${STRAT_TSV}"
    done
fi

# Coverage strata per coverage level
for COV in "${COVERAGES_WGS[@]}"; do
    COV_STRAT_DIR="${STRAT_DIR}/coverage_${COV}x"
    if [[ -d "${COV_STRAT_DIR}" ]]; then
        for cov_bed in "${COV_STRAT_DIR}"/*.bed; do
            name="cov${COV}x_$(basename "${cov_bed}" .bed)"
            echo "${name}	stratification/coverage_${COV}x/$(basename "${cov_bed}")" >> "${STRAT_TSV}"
        done
    fi
done

# ClinVar (will be appended by prepare_clinvar.sh if run)

echo "  Written: ${STRAT_TSV}"
echo "  Entries: $(wc -l < "${STRAT_TSV}")"

# =============================================================================
# 4) GIAB repeat / mappability / segdup strata
# =============================================================================

echo ""
echo "--- Step 4: GIAB repeat/mappability/segdup strata ---"

bash "${SCRIPT_DIR}/prepare_giab_strata.sh"

echo ""
echo "=========================================="
echo " Stratification preparation complete"
echo "=========================================="
echo "  BED files: ${STRAT_DIR}/"
echo "  TSV:       ${STRAT_TSV}"
echo ""
echo "Next: bash evaluation/prepare_clinvar.sh"
