#!/bin/bash
# =============================================================================
# 07_snpeff_annotate.sh
# Chú giải biến thể FN/FP của từng caller bằng SnpEff + SnpSift
# Thêm thông tin: dbSNP, dbNSFP (SIFT, PolyPhen2, PROVEAN, MutationTaster,
#                  CADD, MetaLR), ClinVar
# =============================================================================

set -euo pipefail

# --- Đường dẫn ---
BENCHMARK_DIR="results/benchmarks/rtg_vcfeval"
ANNOT_DIR="results/annotation"

# SnpEff / SnpSift
SNPEFF_JAR="$HOME/snpEff/snpEff.jar"
SNPSIFT_JAR="$HOME/snpEff/SnpSift.jar"
SNPEFF_DB="GRCh38.105"

# Databases
DBNSFP="$HOME/snpEff/data/dbNSFP4.4a.chr22.txt.gz"
CLINVAR_VCF="$HOME/snpEff/data/clinvar/clinvar.vcf.gz"
DBSNP_VCF="data/reference/dbsnp138.hg38.chr22.vcf.gz"

CALLERS="gatk deepvariant strelka2 freebayes"
CATEGORIES="fn fp"

# Các trường dbNSFP cần trích xuất
DBNSFP_FIELDS="SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,PROVEAN_score,PROVEAN_pred,MutationTaster_score,MutationTaster_pred,CADD_phred,MetaLR_score,MetaLR_pred,1000Gp3_AF"

for CALLER in $CALLERS; do
    echo "========================================"
    echo "  Annotating: $CALLER"
    echo "========================================"

    mkdir -p "$ANNOT_DIR/$CALLER"

    for CAT in $CATEGORIES; do
        INPUT="$BENCHMARK_DIR/$CALLER/${CAT}.vcf.gz"
        PREFIX="$ANNOT_DIR/$CALLER/${CAT}"

        echo "--- $CALLER / $CAT ---"

        # Bước 1: SnpEff — chú giải chức năng (missense, nonsense, splice, ...)
        java -Xmx4G -jar "$SNPEFF_JAR" ann \
            -noStats -no-downstream -no-upstream -no-intergenic \
            "$SNPEFF_DB" \
            "$INPUT" \
        | bgzip -c > "${PREFIX}.snpeff.vcf.gz"
        tabix -f -p vcf "${PREFIX}.snpeff.vcf.gz"

        # Bước 2: SnpSift annotate dbSNP — thêm rsID
        java -Xmx4G -jar "$SNPSIFT_JAR" annotate \
            "$DBSNP_VCF" \
            "${PREFIX}.snpeff.vcf.gz" \
        | bgzip -c > "${PREFIX}.dbsnp.vcf.gz"
        tabix -f -p vcf "${PREFIX}.dbsnp.vcf.gz"

        # Bước 3: SnpSift dbnsfp — thêm điểm SIFT, PolyPhen2, PROVEAN, MutationTaster, CADD, MetaLR, 1000Gp3_AF
        java -Xmx4G -jar "$SNPSIFT_JAR" dbnsfp \
            -db "$DBNSFP" \
            -f "$DBNSFP_FIELDS" \
            "${PREFIX}.dbsnp.vcf.gz" \
        | bgzip -c > "${PREFIX}.dbnsfp.vcf.gz"
        tabix -f -p vcf "${PREFIX}.dbnsfp.vcf.gz"

        # Bước 4: SnpSift annotate ClinVar — thêm clinical significance
        java -Xmx4G -jar "$SNPSIFT_JAR" annotate \
            "$CLINVAR_VCF" \
            "${PREFIX}.dbnsfp.vcf.gz" \
        | bgzip -c > "${PREFIX}.annotated.vcf.gz"
        tabix -f -p vcf "${PREFIX}.annotated.vcf.gz"

        # Dọn file trung gian
        rm -f "${PREFIX}.snpeff.vcf.gz"* "${PREFIX}.dbsnp.vcf.gz"* "${PREFIX}.dbnsfp.vcf.gz"*

        echo "  -> Output: ${PREFIX}.annotated.vcf.gz"
    done
done

echo ""
echo "Done. Tất cả file annotated nằm trong $ANNOT_DIR/{caller}/"
