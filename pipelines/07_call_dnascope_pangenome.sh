#!/bin/bash
# Sentieon DNAscope Pangenome — WGS mode from simulated FASTQs.

set -euo pipefail

COV="${1:?Usage: $0 <coverage>  (e.g. 10, 20, 30, 50)}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/_sentieon_common.sh"

CALLER="dnascope_pangenome"
TIMEFILE="${LOG_DIR}/${COV}x/time/dnascope_pangenome.time"

sentieon_set_fastq_paths "${COV}" "wgs"
sentieon_skip_if_no_license "DNAscope Pangenome"
sentieon_require_command docker
sentieon_require_file "${FASTQ_R1}"
sentieon_require_file "${FASTQ_R2}"
sentieon_require_file "${REF_FASTA}"
sentieon_require_file "${DBSNP}"
sentieon_require_file "${DNASCOPE_WGS_MODEL}/bwa.model"
sentieon_require_file "${DNASCOPE_WGS_MODEL}/dnascope.model"
sentieon_require_prefix "${PANGENOME_INDEX}"

sentieon_prepare_layout "${COV_SUFFIX}" "${CALLER}"
sentieon_init_license_args
sentieon_set_pcr_arg

TMPDIR_PG="${OUT_DIR}/tmp"
mkdir -p "${TMPDIR_PG}"

ABS_REF_DIR=$(sentieon_abs_dir "${REF_DIR}")
ABS_SIM_DIR=$(sentieon_abs_dir "${SIM_DIR}")
ABS_OUT_DIR=$(sentieon_abs_dir "${OUT_DIR}")
ABS_TMP_DIR=$(sentieon_abs_dir "${TMPDIR_PG}")
ABS_MODEL_DIR=$(sentieon_abs_dir "${DNASCOPE_WGS_MODEL}")
ABS_PANGENOME_DIR=$(sentieon_abs_dir "$(dirname "${PANGENOME_INDEX}")")

REF_BASENAME=$(basename "${REF_FASTA}")
DBSNP_BASENAME=$(basename "${DBSNP}")
PANGENOME_BASENAME=$(basename "${PANGENOME_INDEX}")
R1_BASENAME=$(basename "${FASTQ_R1}")
R2_BASENAME=$(basename "${FASTQ_R2}")
RG_ID="${SAMPLE_NAME}_${COV}x_pangenome"
RG="@RG\\tID:${RG_ID}\\tSM:${SAMPLE_NAME}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"

echo "=== Running Sentieon DNAscope Pangenome (WGS, supplementary track) ==="

/usr/bin/time -v -o "${TIMEFILE}" \
    docker run --rm \
    "${SENTIEON_DOCKER_LICENSE_ARGS[@]}" \
    -v "${ABS_REF_DIR}:/ref:ro" \
    -v "${ABS_SIM_DIR}:/fastq:ro" \
    -v "${ABS_OUT_DIR}:/output" \
    -v "${ABS_TMP_DIR}:/work" \
    -v "${ABS_MODEL_DIR}:/model:ro" \
    -v "${ABS_PANGENOME_DIR}:/pangenome:ro" \
    "${SENTIEON_IMAGE}" \
    bash -lc "
        set -euo pipefail

        ( sentieon bwa-mem2-pangenome mem \
            -R '${RG}' \
            -t ${THREADS} \
            -K 10000000 \
            -x /model/bwa.model \
            /pangenome/${PANGENOME_BASENAME} \
            /fastq/${R1_BASENAME} /fastq/${R2_BASENAME} ) | \
        sentieon util sort \
            -r /ref/${REF_BASENAME} \
            -o /work/sorted.bam \
            -t ${THREADS} \
            --sam2bam -i -

        sentieon driver \
            -t ${THREADS} \
            -i /work/sorted.bam \
            --algo LocusCollector \
            --fun score_info /work/score.txt

        sentieon driver \
            -t ${THREADS} \
            -i /work/sorted.bam \
            --algo Dedup \
            --score_info /work/score.txt \
            --metrics /output/${PREFIX}_PANGENOME_dedup_metrics.txt \
            /work/deduped.bam

        sentieon driver \
            -r /ref/${REF_BASENAME} \
            -t ${THREADS} \
            -i /work/deduped.bam \
            --algo DNAscope ${SENTIEON_PCR_ARG} \
            --model /model/dnascope.model \
            -d /ref/${DBSNP_BASENAME} \
            /output/${PREFIX}_DNASCOPE_PG_TMP.vcf.gz

        sentieon driver \
            -r /ref/${REF_BASENAME} \
            -t ${THREADS} \
            --algo DNAModelApply \
            --model /model/dnascope.model \
            -v /output/${PREFIX}_DNASCOPE_PG_TMP.vcf.gz \
            /output/${PREFIX}_DNASCOPE_PANGENOME.vcf.gz
    "

sentieon_log_metrics "dnascope_pangenome" "DNAscope_Pangenome" "${TIMEFILE}"

rm -rf "${TMPDIR_PG}"
rm -f "${OUT_DIR}/${PREFIX}_DNASCOPE_PG_TMP.vcf.gz" "${OUT_DIR}/${PREFIX}_DNASCOPE_PG_TMP.vcf.gz.tbi"

echo
echo "DNAscope Pangenome done: ${OUT_DIR}/${PREFIX}_DNASCOPE_PANGENOME.vcf.gz"
echo "Metrics: ${METRICS}"
