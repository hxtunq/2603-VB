#!/usr/bin/env bash
# =============================================================================
# Evaluate variant callers against truth set using RTG vcfeval (directly).
# Replaces hap.py-based evaluation — no Docker dependency.
#
# For each coverage × caller, produces:
#   results/eval/vcfeval/{COV}x/{caller}/  →  fn.vcf.gz, fp.vcf.gz, tp.vcf.gz,
#                                              tp-baseline.vcf.gz, summary.txt
#
# Aggregated stats written to: results/eval/vcfeval_all_stats.tsv
#
# Usage:
#   bash evaluation/eval_vcfeval.sh                  # all callers × all coverages
#   WGS_CALLERS="gatk deepvariant" bash evaluation/eval_vcfeval.sh   # subset
# =============================================================================

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "${SCRIPT_DIR}/../config/config.sh"

# ---------- Callers ----------
WGS_CALLERS="${WGS_CALLERS:-gatk deepvariant strelka2 freebayes dnascope}"

# ---------- VCF naming (must match pipelines/*.sh) ----------
get_vcf() {
    local caller="$1" cov="$2"
    local variant_dir="${VARIANT_DIR}/${cov}x"
    case "${caller}" in
        gatk)            echo "${variant_dir}/gatk/SS_${CHR_TO_USE}_HC_${cov}x_WGS.vcf.gz" ;;
        deepvariant)     echo "${variant_dir}/deepvariant/SS_${CHR_TO_USE}_DV_${cov}x_WGS.vcf.gz" ;;
        strelka2)        echo "${variant_dir}/strelka2/SS_${CHR_TO_USE}_STRELKA_${cov}x_WGS.vcf.gz" ;;
        freebayes)       echo "${variant_dir}/freebayes/SS_${CHR_TO_USE}_FB_${cov}x_WGS.vcf.gz" ;;
        dnascope)        echo "${variant_dir}/dnascope/SS_${CHR_TO_USE}_DNASCOPE_${cov}x_WGS.vcf.gz" ;;
        *)               echo "ERROR: unknown caller: ${caller}" >&2; return 1 ;;
    esac
}

caller_alias() {
    case "$1" in
        gatk)        echo "HC" ;;
        deepvariant) echo "DV" ;;
        strelka2)    echo "ST" ;;
        freebayes)   echo "FB" ;;
        dnascope)    echo "DS" ;;
        *)           echo "$1" ;;
    esac
}

# ---------- Prepare truth VCF with GT column ----------
prepare_truth() {
    local truth_prepared="${EVAL_DIR}/vcfeval/truth_prepared.vcf.gz"

    if [[ -f "${truth_prepared}" ]]; then
        echo "[TRUTH] Already prepared: ${truth_prepared}"
        return
    fi

    mkdir -p "$(dirname "${truth_prepared}")"
    echo "[TRUTH] Preparing truth VCF with GT column..."

    # simuG produces sites-only VCFs; add ##FORMAT=GT header + GT:1/1
    zcat "${TRUTH_VCF}" | awk 'BEGIN{OFS="\t"; ff=0}
      /^##fileformat=/ {ff=1; print; next}
      /^##/ {print; next}
      /^#CHROM/ {
        if(ff==0) print "##fileformat=VCFv4.2"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print $0,"FORMAT","TRUTH"
        next
      }
      {print $0,"GT","1/1"}
    ' | bcftools norm -f "${REF_FASTA}" -m -both \
      | bgzip -c > "${truth_prepared}"

    tabix -f -p vcf "${truth_prepared}"
    echo "[TRUTH] Done: ${truth_prepared}"
}

# ---------- PASS-filter + normalize a caller VCF ----------
prepare_query() {
    local vcf="$1"
    local norm="${vcf%.vcf.gz}.pass_norm.vcf.gz"

    if [[ -f "${norm}" ]]; then
        echo "${norm}"
        return
    fi

    if [[ ! -f "${vcf}" ]]; then
        echo ""
        return
    fi

    echo "  [PREP] PASS-filter + normalize: $(basename "${vcf}")" >&2
    bcftools view -f 'PASS,.' "${vcf}" \
        | bcftools norm -m -both -f "${REF_FASTA}" \
        | bgzip -c > "${norm}"
    tabix -f -p vcf "${norm}"
    echo "${norm}"
}

# ---------- Create SDF template if needed ----------
ensure_sdf() {
    if [[ -d "${RTG_SDF}" ]]; then
        return
    fi
    echo "[SDF] Creating SDF template: ${RTG_SDF}"
    rtg format -o "${RTG_SDF}" "${REF_FASTA}"
}

# ---------- Run vcfeval for one caller × one coverage ----------
run_vcfeval_one() {
    local caller="$1" cov="$2"
    local alias
    alias=$(caller_alias "${caller}")

    local vcf_raw
    vcf_raw=$(get_vcf "${caller}" "${cov}")

    local vcf_norm
    vcf_norm=$(prepare_query "${vcf_raw}")

    if [[ -z "${vcf_norm}" ]]; then
        echo "  SKIP: ${alias} @ ${cov}x — VCF not found: ${vcf_raw}"
        return
    fi

    local outdir="${EVAL_DIR}/vcfeval/${cov}x/${caller}"

    if [[ -f "${outdir}/summary.txt" ]]; then
        echo "  EXISTS: ${alias} @ ${cov}x — skipping"
        return
    fi

    # Remove partial output if present
    [[ -d "${outdir}" ]] && rm -rf "${outdir}"

    echo "  Running: ${alias} @ ${cov}x"
    echo "    Query: ${vcf_norm}"

    rtg vcfeval \
        --baseline "${EVAL_DIR}/vcfeval/truth_prepared.vcf.gz" \
        --calls "${vcf_norm}" \
        --template "${RTG_SDF}" \
        --bed-regions "${TRUTH_BED}" \
        --output "${outdir}" \
        --threads "${THREADS}"

    echo "  Done:  ${outdir}/summary.txt"
}

# ---------- Aggregate summary.txt → TSV ----------
aggregate_stats() {
    local tsv="$1"
    local vcfeval_root="${EVAL_DIR}/vcfeval"

    echo "Coverage	Caller	CallerAlias	Threshold	True-pos-baseline	True-pos-call	False-pos	False-neg	Precision	Sensitivity	F-measure" > "${tsv}"

    find "${vcfeval_root}" -name "summary.txt" | sort | while IFS= read -r summary; do
        # path: .../vcfeval/{COV}x/{caller}/summary.txt
        local rel="${summary#${vcfeval_root}/}"
        local cov_dir="${rel%%/*}"        # e.g. "30x"
        local caller; caller=$(echo "${rel}" | cut -d'/' -f2)
        local alias; alias=$(caller_alias "${caller}")

        # Read the data lines (skip header lines starting with spaces or dashes)
        tail -n +3 "${summary}" | while IFS= read -r line; do
            # summary.txt format:  Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
            # Each line may be indented; trim leading spaces
            local trimmed; trimmed=$(echo "${line}" | sed 's/^[[:space:]]*//')
            [[ -z "${trimmed}" ]] && continue
            [[ "${trimmed}" == -* ]] && continue

            printf '%s\t%s\t%s\t%s\n' "${cov_dir}" "${caller}" "${alias}" "$(echo "${trimmed}" | tr -s ' ' '\t')"
        done
    done >> "${tsv}"

    echo "Stats aggregated: ${tsv}"
    echo "Total lines: $(wc -l < "${tsv}")"
}

# =============================================================================
# MAIN
# =============================================================================

echo "=========================================="
echo " RTG VCFeval Evaluation"
echo "=========================================="
echo "Callers:   ${WGS_CALLERS}"
echo "Coverages: ${COVERAGES_WGS[*]}"
echo ""

# Pre-flight checks
command -v rtg &>/dev/null || { echo "ERROR: 'rtg' not found in PATH"; exit 1; }
command -v bcftools &>/dev/null || { echo "ERROR: 'bcftools' not found in PATH"; exit 1; }
command -v bgzip &>/dev/null || { echo "ERROR: 'bgzip' not found in PATH"; exit 1; }
command -v tabix &>/dev/null || { echo "ERROR: 'tabix' not found in PATH"; exit 1; }
[[ -f "${TRUTH_VCF}" ]] || { echo "ERROR: Truth VCF not found: ${TRUTH_VCF}"; exit 1; }
[[ -f "${REF_FASTA}" ]] || { echo "ERROR: Reference FASTA not found: ${REF_FASTA}"; exit 1; }
[[ -f "${TRUTH_BED}" ]] || { echo "ERROR: BED not found: ${TRUTH_BED}"; exit 1; }

ensure_sdf
prepare_truth

# Run vcfeval for each coverage × caller
for COV in "${COVERAGES_WGS[@]}"; do
    echo ""
    echo "--- Coverage: ${COV}x ---"
    for caller in ${WGS_CALLERS}; do
        run_vcfeval_one "${caller}" "${COV}"
    done
done

# Aggregate all results
echo ""
echo "--- Aggregating stats ---"
STATS_TSV="${EVAL_DIR}/vcfeval_all_stats.tsv"
aggregate_stats "${STATS_TSV}"

echo ""
echo "=========================================="
echo " Evaluation complete"
echo "=========================================="
echo "Per-caller results: ${EVAL_DIR}/vcfeval/{COV}x/{caller}/"
echo "Aggregated stats:   ${STATS_TSV}"
echo ""
echo "Next steps:"
echo "  python visualization/plot_summary.py ${STATS_TSV} results/plots"
echo "  python visualization/plot_coverage_comparison.py ${STATS_TSV} results/plots"
