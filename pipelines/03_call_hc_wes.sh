#!/bin/bash
# GATK HaplotypeCaller + Hard Filtering — WES Mode
# Same as 03_call_hc.sh but with --intervals for exome targets
# CNN/NVScore removed — GATK 4.6 scorevariants package unavailable

set -e

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50, 100, 200)}"

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/${COV}x/bam_path.sh"

CALLER="gatk_wes"
OUT_DIR="${VARIANT_DIR}/${COV}x/${CALLER}"
TIMEDIR="${LOG_DIR}/${COV}x/time"
mkdir -p "${OUT_DIR}" "${TIMEDIR}"

METRICS="${LOG_DIR}/${COV}x/benchmark_metrics.tsv"
if [[ ! -f "${METRICS}" ]]; then
    echo -e "Caller\tPipeline\tWallClock_sec\tCPU_percent\tMaxRSS_kB" > "${METRICS}"
fi

log_metrics() {
    local caller="$1" pipeline="$2" timefile="$3"
    local wall cpu rss
    wall=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "$timefile")
    cpu=$(awk -F': ' '/Percent of CPU this job got/{print $2}' "$timefile" | tr -d '%')
    rss=$(awk -F': ' '/Maximum resident set size \(kbytes\)/{print $2}' "$timefile")
    local secs
    secs=$(echo "$wall" | awk -F: '{if(NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}')
    echo -e "${caller}\t${pipeline}\t${secs}\t${cpu}\t${rss}" >> "${METRICS}"
    echo "  [METRICS] ${pipeline}: ${wall} wall, ${cpu}% CPU, MaxRSS=${rss} kB"
}

# Validate exome BED
if [[ ! -f "${EXOME_BED}" ]]; then
    echo "ERROR: EXOME_BED not found: ${EXOME_BED}"
    exit 1
fi

# --- BQSR (not timed) ---
echo "Running BaseRecalibrator (WES)..."
gatk --java-options "${JAVA_OPTS}" BaseRecalibrator \
    -R "${REF_FASTA}" \
    -I "${FINAL_BAM}" \
    --known-sites "${DBSNP}" \
    --known-sites "${KNOWN_INDELS}" \
    -L "${EXOME_BED}" \
    -O "${OUT_DIR}/${PREFIX}.recal.table"

gatk --java-options "${JAVA_OPTS}" ApplyBQSR \
    -R "${REF_FASTA}" \
    -I "${FINAL_BAM}" \
    -bqsr "${OUT_DIR}/${PREFIX}.recal.table" \
    -O "${OUT_DIR}/${PREFIX}.recal.bam"

# ==========================================================================
# HaplotypeCaller (WES mode with intervals)
# ==========================================================================
echo "=== HaplotypeCaller (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/hc_wes_calling.time" \
    gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
    -R "${REF_FASTA}" \
    -I "${OUT_DIR}/${PREFIX}.recal.bam" \
    -O "${OUT_DIR}/${PREFIX}_HC.RAW.vcf" \
    -L "${EXOME_BED}" \
    --standard-min-confidence-threshold-for-calling "${GATK_STAND_CALL_CONF}" \
    --native-pair-hmm-threads "${THREADS}"

log_metrics "gatk_wes" "HaplotypeCaller_WES" "${TIMEDIR}/hc_wes_calling.time"

# ==========================================================================
# Hard Filter (GATK Best Practices)
# ==========================================================================
echo "=== Filter: Hard Filter (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/hc_wes_filter_hard.time" bash -c "
    set -e
    gatk SelectVariants -V '${OUT_DIR}/${PREFIX}_HC.RAW.vcf' \
        -select-type SNP -select-type MIXED \
        -O '${OUT_DIR}/${PREFIX}.HC.snps.vcf'

    gatk VariantFiltration -V '${OUT_DIR}/${PREFIX}.HC.snps.vcf' \
        -filter 'QD < 2.0' --filter-name 'QD2' \
        -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
        -filter 'SOR > 3.0' --filter-name 'SOR3' \
        -filter 'FS > 60.0' --filter-name 'FS60' \
        -filter 'MQ < 40.0' --filter-name 'MQ40' \
        -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
        -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \
        -O '${OUT_DIR}/${PREFIX}.HC.snps.flt.vcf'

    gatk SelectVariants -V '${OUT_DIR}/${PREFIX}_HC.RAW.vcf' \
        -select-type INDEL -select-type MIXED \
        -O '${OUT_DIR}/${PREFIX}.HC.indels.vcf'

    gatk VariantFiltration -V '${OUT_DIR}/${PREFIX}.HC.indels.vcf' \
        -filter 'QD < 2.0' --filter-name 'QD2' \
        -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
        -filter 'FS > 200.0' --filter-name 'FS200' \
        -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \
        -O '${OUT_DIR}/${PREFIX}.HC.indels.flt.vcf'

    gatk MergeVcfs \
        -I '${OUT_DIR}/${PREFIX}.HC.snps.flt.vcf' \
        -I '${OUT_DIR}/${PREFIX}.HC.indels.flt.vcf' \
        -O '${OUT_DIR}/${PREFIX}_HC_HARDFILTER.vcf'
"
log_metrics "gatk_wes" "Filter_HardFilter_WES" "${TIMEDIR}/hc_wes_filter_hard.time"

echo ""
echo "HC WES done. Metrics in ${METRICS}:"
echo "  HaplotypeCaller_WES   = calling with -L ${EXOME_BED}"
echo "  Filter_HardFilter_WES → ${PREFIX}_HC_HARDFILTER.vcf"
