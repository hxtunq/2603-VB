#!/bin/bash
# GATK HaplotypeCaller + CNN Scoring + Hard Filtering — WES Mode
# Same as 03_call_hc.sh but with --intervals for exome targets
# Uses EXOME_BED to restrict calling to exome regions

set -e

source "$(dirname "$0")/../config/config.sh"
source "${PREPROC_DIR}/bam_path.sh"

CALLER="gatk_wes"
OUT_DIR="${VARIANT_DIR}/${CALLER}"
TIMEDIR="${LOG_DIR}/time"
mkdir -p "${OUT_DIR}" "${TIMEDIR}"

METRICS="${LOG_DIR}/benchmark_metrics.tsv"
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
# HaplotypeCaller (measured separately — WES mode with intervals)
# ==========================================================================
echo "=== HaplotypeCaller (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/hc_wes_calling.time" \
    gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
    -R "${REF_FASTA}" \
    -I "${OUT_DIR}/${PREFIX}.recal.bam" \
    -bamout "${OUT_DIR}/${PREFIX}.HC.bam" \
    -O "${OUT_DIR}/${PREFIX}_HC.RAW.vcf" \
    -L "${EXOME_BED}" \
    --standard-min-confidence-threshold-for-calling "${GATK_STAND_CALL_CONF}" \
    --native-pair-hmm-threads "${THREADS}"

log_metrics "gatk_wes" "HaplotypeCaller_WES" "${TIMEDIR}/hc_wes_calling.time"

# ==========================================================================
# Filtering 1: CNN 1D
# ==========================================================================
echo "=== Filter: CNN 1D (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/hc_wes_filter_cnn1d.time" bash -c "
    set -e
    gatk --java-options '${JAVA_OPTS}' CNNScoreVariants \
        -R '${REF_FASTA}' \
        -V '${OUT_DIR}/${PREFIX}_HC.RAW.vcf' \
        -O '${OUT_DIR}/${PREFIX}.HC.CNN1.vcf' \
        -tensor-type reference

    gatk --java-options '${JAVA_OPTS}' FilterVariantTranches \
        -V '${OUT_DIR}/${PREFIX}.HC.CNN1.vcf' \
        --output '${OUT_DIR}/${PREFIX}_HC_1DCNN.vcf' \
        --info-key CNN_1D \
        --snp-tranche ${CNN_SNP_TRANCHE} --indel-tranche ${CNN_INDEL_TRANCHE} \
        --resource '${KNOWN_INDELS}' \
        --resource '${KNOWN_SNPS}' \
        --resource '${DBSNP}'
"
log_metrics "gatk_wes" "Filter_CNN_1D_WES" "${TIMEDIR}/hc_wes_filter_cnn1d.time"

# ==========================================================================
# Filtering 2: CNN 2D
# ==========================================================================
echo "=== Filter: CNN 2D (WES) ==="
/usr/bin/time -v -o "${TIMEDIR}/hc_wes_filter_cnn2d.time" bash -c "
    set -e
    gatk --java-options '${JAVA_OPTS}' CNNScoreVariants \
        -I '${OUT_DIR}/${PREFIX}.HC.bam' \
        -R '${REF_FASTA}' \
        -V '${OUT_DIR}/${PREFIX}.HC.CNN1.vcf' \
        -O '${OUT_DIR}/${PREFIX}.HC.CNN2.vcf' \
        -tensor-type read_tensor

    gatk --java-options '${JAVA_OPTS}' FilterVariantTranches \
        -V '${OUT_DIR}/${PREFIX}.HC.CNN2.vcf' \
        --output '${OUT_DIR}/${PREFIX}_HC_2DCNN.vcf' \
        --info-key CNN_2D \
        --snp-tranche ${CNN_SNP_TRANCHE} --indel-tranche ${CNN_INDEL_TRANCHE} \
        --resource '${KNOWN_INDELS}' \
        --resource '${KNOWN_SNPS}' \
        --resource '${DBSNP}'
"
log_metrics "gatk_wes" "Filter_CNN_2D_WES" "${TIMEDIR}/hc_wes_filter_cnn2d.time"

# ==========================================================================
# Filtering 3: Hard Filter
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
echo "  HaplotypeCaller_WES  = calling with -L ${EXOME_BED}"
echo "  Filter_CNN_1D_WES    → ${PREFIX}_HC_1DCNN.vcf"
echo "  Filter_CNN_2D_WES    → ${PREFIX}_HC_2DCNN.vcf"
echo "  Filter_HardFilter_WES → ${PREFIX}_HC_HARDFILTER.vcf"
