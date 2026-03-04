## Quy trình thực hiện

### Phần I. Tạo cấu trúc folder

```bash
git clone https://github.com/bitschif/variant-benchmarking.git
cd variant-benchmarking/
```

```bash
# pwd: variant-benchmarking

mkdir -p data/reference
mkdir -p data/simulated
mkdir -p results/preprocessing
mkdir -p results/variants/{gatk,deepvariant,strelka2,freebayes,dnascope}
mkdir -p results/variants/{gatk_wes,deepvariant_wes,strelka2_wes,freebayes_wes,dnascope_wes}
mkdir -p results/eval
mkdir -p results/plots
mkdir -p logs/time
```

### Phần II. Download và chuẩn bị dữ liệu

#### 2.1. Download Reference Genome (chr22 - hg38)

```bash
#pwd: variant-benchmarking/data/reference

# download từ UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# giải nén
gunzip chr22.fa.gz

# index reference
samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

Output:
- `data/reference/chr22.fa` — reference genome
- `data/reference/chr22.fa.fai` — samtools index
- `data/reference/chr22.fa.bwt` (+ `.amb`, `.ann`, `.pac`, `.sa`) — BWA index
- `data/reference/chr22.dict` — sequence dictionary

#### 2.2. Download Known Sites cho BQSR

Known sites giúp cải thiện chất lượng Base Quality Score Recalibration.

```bash
# pwd: variant-benchmarking/data/reference

# === dbSNP ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Extract chr22 và đảm bảo chromosome naming là "chr22"
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.vcf.gz

bcftools view -r chr22 -Oz -o dbsnp138.hg38.chr22.vcf.gz dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.chr22.vcf.gz

# === Mills and 1000G Gold Standard Indels ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# extract chr22
bcftools view -r chr22 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

# === 1000G Phase 1 SNPs ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# extract chr22
bcftools view -r chr22 1000G_phase1.snps.high_confidence.hg38.vcf.gz -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz

# xoá file không sử dụng đến
rm -f Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.dbsnp138.vcf.idx
rm -f Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
rm -f 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

**Một số lưu ý:**

```bash
# pwd: variant-benchmarking/data/reference

# reference genome từ UCSC sử dụng format "chr22"
# tất cả VCF files phải có tên chromosome khớp với tên ở file reference

# lệnh kiểm tra tên chromosome trong các file VCF
bcftools view -h file.vcf.gz | grep "^##contig"
tabix -l gile.vcf.gz # ở đồ án và repository này lấy tên là "chr22" làm chuẩn

# nếu gặp lỗi "chromosome not found", kiểm tra và đổi lại tên chromosome
# cách đổi tên "22" thành "chr22" nếu cần
bcftools annotate --rename-chrs chr_map.txt input.vcf.gz -Oz -o output.vcf.gz

# tạo file BED chứa vùng non-N của chromosome
python3 << 'EOF'
from Bio import SeqIO
import re

fasta = "chr22.fa"
output = "chr22_non_N_regions.bed"

with open(output, 'w') as out:
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq).upper()
        chrom = record.id
        for m in re.finditer(r'[ATGC]+', seq):
            start = m.start()
            end = m.end()
            if end - start >= 1000:
                out.write(f"{chrom}\t{start}\t{end}\n")

print("Done")
EOF
```

Output:
- `data/reference/dbsnp138.hg38.chr22.vcf.gz` (+`.tbi`) — dbSNP chr22
- `data/reference/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz` (+`.tbi`) — Mills indels chr22
- `data/reference/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz` (+`.tbi`) — 1000G SNPs chr22
- `data/reference/chr22_non_N_regions.bed` — vùng non-N cho VCFeval

#### 2.3. Tạo đột biến SNPs và INDELs với simuG

Thông số được chỉnh theo dữ liệu thực tế (chr22 ~50 Mb, ~20K variants cho GIAB samples).

> **Lưu ý:** KHÔNG dùng `-titv_ratio` để tránh bias benchmark. Để simuG tự random
> generate SNPs với tỉ lệ transition/transversion tự nhiên.

```bash
# pwd: variant-benchmarking/data

# clone simuG về folder <data>
git clone https://github.com/yjx1217/simuG.git

# dùng simuG tạo đột biến — thông số realistic
perl simuG/simuG.pl \
  -refseq reference/chr22.fa \
  -snp_count 18000 \
  -indel_count 2500 \
  -indel_min_len 1 \
  -indel_max_len 30 \
  -prefix simulated/SIMULATED_SAMPLE_chr22
```

Output:
- `data/simulated/SIMULATED_SAMPLE_chr22.simseq.genome.fa` — genome đã chứa đột biến
- `data/simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.vcf` — danh sách SNP đã tạo (~18K)
- `data/simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.vcf` — danh sách INDEL đã tạo (~2.5K)

#### 2.4. Merge các file VCF output của simuG tạo truth VCF

```bash
# pwd: variant-benchmarking/data

# đảm bảo có .fai từ reference
samtools faidx reference/chr22.fa

# reheader từng file VCF với contig info
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.vcf

bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.vcf

# concat sau khi đã có header chuẩn
bcftools concat \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf | \
bcftools sort -Oz -o simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

tabix -p vcf simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz
```

Output: `data/simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz` (+`.tbi`) — truth VCF chứa toàn bộ SNP + INDEL

#### 2.5. Tạo file fastq giả lập bằng ART Illumina — Multi-coverage

Tạo nhiều mức coverage khác nhau từ CÙNG MỘT mutated genome:
- **WGS coverage:** 10x, 20x, 30x, 50x (tương tự GIAB WGS 22–37x)
- **WES coverage:** 100x, 200x (tương tự GIAB WES 183–249x, đánh giá trên exome BED)

```bash
# pwd: variant-benchmarking/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

# Multi-coverage loop
for COV in 10 20 30 50 100 200; do
  echo "=== Generating ${COV}x coverage ==="

  art_illumina \
      -ss HS25 \
      -i "${MUTATED_FASTA}" \
      -p \
      -l 150 \
      -f ${COV} \
      -m 350 \
      -s 50 \
      -rs 42 \
      -o "${SIM_DIR}/${PREFIX}_${COV}x_" \
      -na

  # đổi tên và nén
  mv "${SIM_DIR}/${PREFIX}_${COV}x_1.fq" "${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq"
  mv "${SIM_DIR}/${PREFIX}_${COV}x_2.fq" "${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq"

  echo "Done: ${SIM_DIR}/${PREFIX}_${COV}x_R{1,2}.fastq.gz"
done
```

Output (per coverage level):
- `data/simulated/SIMULATED_SAMPLE_chr22_{COV}x_R1.fastq.gz` — paired-end read 1
- `data/simulated/SIMULATED_SAMPLE_chr22_{COV}x_R2.fastq.gz` — paired-end read 2

#### 2.6. Download BED files cho exome targets và stratification

Lấy từ AstraZeneca-NGS reference_data (hg38), extract chr22.

```bash
# pwd: variant-benchmarking/data/reference

# === CDS regions (canonical protein-coding) ===
wget -O CDS-canonical.bed \
  https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/CDS-canonical.bed

# Extract chr22 only
grep -w "chr22" CDS-canonical.bed > CDS-canonical.chr22.bed
rm CDS-canonical.bed

# === Exome target BED (Agilent SureSelect V6) ===
wget -O Exome-Agilent_V6.bed \
  https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/Exome-Agilent_V6.bed

# Extract chr22 only
grep -w "chr22" Exome-Agilent_V6.bed > Exome-Agilent_V6.chr22.bed
rm Exome-Agilent_V6.bed

# === Mappability tricky regions ===
wget -O umap_k100_mappability.bed.gz \
  https://github.com/AstraZeneca-NGS/reference_data/raw/master/hg38/tricky_regions/umap_k100_mappability.bed.gz
wget -O umap_k100_mappability.bed.gz.tbi \
  https://github.com/AstraZeneca-NGS/reference_data/raw/master/hg38/tricky_regions/umap_k100_mappability.bed.gz.tbi

# Extract chr22 from mappability
tabix umap_k100_mappability.bed.gz chr22 | bgzip > umap_k100_mappability.chr22.bed.gz
tabix -p bed umap_k100_mappability.chr22.bed.gz
rm umap_k100_mappability.bed.gz umap_k100_mappability.bed.gz.tbi

# === Tạo stratification TSV cho hap.py ===
cat > stratification_chr22.tsv <<'EOF'
CDS	CDS-canonical.chr22.bed
Exome_Targets	Exome-Agilent_V6.chr22.bed
Low_Mappability	umap_k100_mappability.chr22.bed.gz
EOF
```

Output:
- `data/reference/CDS-canonical.chr22.bed` — CDS regions chr22
- `data/reference/Exome-Agilent_V6.chr22.bed` — exome targets chr22
- `data/reference/umap_k100_mappability.chr22.bed.gz` — low-mappability regions chr22
- `data/reference/stratification_chr22.tsv` — stratification file cho hap.py

### Phần III. Tiền xử lý dữ liệu cho công đoạn gọi biến thể — Multi-coverage

Pipeline chạy lặp cho MỖI mức coverage: FastQC → BWA-MEM → MarkDuplicates → BQSR → Filter

```bash
# pwd: variant-benchmarking

PREFIX="SIMULATED_SAMPLE_chr22"
REF="data/reference/chr22.fa"

for COV in 10 20 30 50 100 200; do
  echo "============================================"
  echo "Processing ${COV}x coverage"
  echo "============================================"

  R1="data/simulated/${PREFIX}_${COV}x_R1.fastq.gz"
  R2="data/simulated/${PREFIX}_${COV}x_R2.fastq.gz"
  OUTDIR="results/preprocessing/${COV}x"
  mkdir -p "${OUTDIR}" "logs/${COV}x"

  # --- 3.1 FastQC ---
  mkdir -p "${OUTDIR}/fastqc"
  fastqc -t 4 -o "${OUTDIR}/fastqc" "${R1}" "${R2}" 2>&1 | tee "logs/${COV}x/fastqc.log"

  # --- 3.2 BWA-MEM ---
  bwa mem -t 4 -M \
    -R "@RG\tID:${PREFIX}_${COV}x\tSM:${PREFIX}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${R1}" "${R2}" \
    2> "logs/${COV}x/bwa_mem.log" | \
    samtools sort -@ 4 -m 2G -o "${OUTDIR}/${PREFIX}_aligned.bam" -
  samtools index "${OUTDIR}/${PREFIX}_aligned.bam"

  # --- 3.3 MarkDuplicates ---
  gatk MarkDuplicates \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_marked.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true 2>&1 | tee "logs/${COV}x/markduplicates.log"

  # --- 3.4 BQSR ---
  gatk BaseRecalibrator \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -R "${REF}" \
    -I "${OUTDIR}/${PREFIX}_marked.bam" \
    --known-sites data/reference/dbsnp138.hg38.chr22.vcf.gz \
    --known-sites data/reference/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
    --known-sites data/reference/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz \
    -O "${OUTDIR}/${PREFIX}_recal.table" 2>&1 | tee "logs/${COV}x/bqsr.log"

  gatk ApplyBQSR \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -R "${REF}" -I "${OUTDIR}/${PREFIX}_marked.bam" \
    --bqsr-recal-file "${OUTDIR}/${PREFIX}_recal.table" \
    -O "${OUTDIR}/${PREFIX}_recal.bam" 2>&1 | tee "logs/${COV}x/applybqsr.log"
  samtools index "${OUTDIR}/${PREFIX}_recal.bam"

  # --- 3.5 Filter by mapping quality ---
  samtools view -@ 4 -b -q 20 -F 1796 "${OUTDIR}/${PREFIX}_recal.bam" | \
    samtools sort -@ 4 -o "${OUTDIR}/${PREFIX}_final.bam" -
  samtools index "${OUTDIR}/${PREFIX}_final.bam"

  # --- 3.6 Coverage stats ---
  samtools stats "${OUTDIR}/${PREFIX}_final.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
  samtools flagstat "${OUTDIR}/${PREFIX}_final.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
  mosdepth -t 4 --by 1000 "${OUTDIR}/${PREFIX}_coverage" "${OUTDIR}/${PREFIX}_final.bam"

  # --- 3.7 Export BAM path ---
  echo "FINAL_BAM=${OUTDIR}/${PREFIX}_final.bam" > "${OUTDIR}/bam_path.sh"

  echo "Done: ${OUTDIR}/${PREFIX}_final.bam"
done
```

Output (per coverage level `{COV}x`):
- `results/preprocessing/{COV}x/SIMULATED_SAMPLE_chr22_final.bam` — BAM cuối cùng
- `results/preprocessing/{COV}x/bam_path.sh` — export path cho downstream scripts

### Phần IV. Gọi biến thể — Multi-coverage, WGS + WES

Chạy tất cả caller cho mỗi mức coverage. WGS callers cho mọi coverage,
WES callers cho coverage phù hợp (50x, 100x, 200x).

#### 4.1 WGS Variant Calling

```bash
# pwd: variant-benchmarking

for COV in 10 20 30 50 100 200; do
  echo "=== Variant Calling: ${COV}x WGS ==="
  export PREPROC_DIR="results/preprocessing/${COV}x"
  export VARIANT_DIR="results/variants/${COV}x"
  export LOG_DIR="logs/${COV}x"
  mkdir -p "${VARIANT_DIR}"/{gatk,deepvariant,strelka2,freebayes,dnascope}
  mkdir -p "${LOG_DIR}/time"

  source "${PREPROC_DIR}/bam_path.sh"

  bash pipelines/03_call_hc.sh                   # GATK HC (3 filters)
  bash pipelines/04_call_dv.sh                    # DeepVariant
  bash pipelines/05_call_strelka.sh               # Strelka2 (with Manta)
  bash pipelines/06_call_freebayes.sh             # FreeBayes
  bash pipelines/07_call_dnascope.sh              # DNAscope (optional)
done
```

#### 4.2 WES Variant Calling (coverage ≥ 50x)

```bash
# pwd: variant-benchmarking

for COV in 50 100 200; do
  echo "=== Variant Calling: ${COV}x WES ==="
  export PREPROC_DIR="results/preprocessing/${COV}x"
  export VARIANT_DIR="results/variants/${COV}x"
  export LOG_DIR="logs/${COV}x"
  mkdir -p "${VARIANT_DIR}"/{gatk_wes,deepvariant_wes,strelka2_wes,freebayes_wes,dnascope_wes}

  source "${PREPROC_DIR}/bam_path.sh"

  bash pipelines/03_call_hc_wes.sh               # GATK HC WES
  bash pipelines/04_call_dv_wes.sh                # DeepVariant WES
  bash pipelines/05_call_strelka_wes.sh           # Strelka2 WES (with Manta --exome)
  bash pipelines/06_call_freebayes_wes.sh         # FreeBayes WES
  bash pipelines/07_call_dnascope_wes.sh          # DNAscope WES (optional)
done
```

Output per coverage:
- `results/variants/{COV}x/gatk/` — GATK HC (3 filter strategies)
- `results/variants/{COV}x/deepvariant/` — DeepVariant
- `results/variants/{COV}x/strelka2/` — Strelka2 (with Manta indel candidates)
- `results/variants/{COV}x/freebayes/` — FreeBayes
- `results/variants/{COV}x/dnascope/` — DNAscope (optional, needs Sentieon license)
- `results/variants/{COV}x/*_wes/` — WES mode scripts (same but with exome targets)

> **Callers:** 5 callers × (WGS + WES) = 10 pipelines per coverage level
> **Total VCFs:** 6 coverages × WGS (5) + 3 coverages × WES (5) = 45+ VCFs

### Phần V. So sánh các công cụ gọi biến thể bằng RTG Tools VCFEval 

#### 5.1 Chạy VCFeval cho 4 file VCF output của 4 công cụ gọi biến thể sau khi gọi biến thể

```bash
# pwd: variant-benchmarking

# khai báo đường dẫn
BED=data/reference/chr22_non_N_regions.bed
REF=data/reference/chr22.fa
TRUTH_RAW=$(ls -1 data/simulated/*_truth.vcf.gz | head -n 1)

mkdir -p results/benchmarks/truth
rm -f results/benchmarks/truth/truth.gt.vcf.gz* results/benchmarks/truth/truth.gt.norm.vcf.gz*

# Tạo truth có FORMAT+GT hợp lệ (có khai báo ##FORMAT cho GT)
zcat "$TRUTH_RAW" | awk 'BEGIN{OFS="\t"; ff=0}
  /^##fileformat=/ {ff=1; print; next}
  /^##/ {print; next}
  /^#CHROM/ {
    if(ff==0) print "##fileformat=VCFv4.2"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print $0,"FORMAT","TRUTH"
    next
  }
  {print $0,"GT","1/1"}
' | bgzip -c > results/benchmarks/truth/truth.gt.vcf.gz

tabix -f -p vcf results/benchmarks/truth/truth.gt.vcf.gz

# Normalize lại truth (có GT) để dùng cho vcfeval
bcftools norm -f "$REF" -m -both results/benchmarks/truth/truth.gt.vcf.gz -Oz \
  -o results/benchmarks/truth/truth.gt.norm.vcf.gz

tabix -f -p vcf results/benchmarks/truth/truth.gt.norm.vcf.gz

# tạo SDF template cho RTG
rtg format -o results/benchmarks/ref/chr22.sdf "$REF"
# check
bcftools query -l results/benchmarks/truth/truth.gt.norm.vcf.gz
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/gatk/*_gatk_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/gatk \
  --threads 4
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/deepvariant/*_deepvariant_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/deepvariant \
  --threads 4
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/strelka2/*_strelka2_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/strelka2 \
  --threads 4
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/freebayes/*_freebayes_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/freebayes \
  --threads 4
```

Output per caller (`results/benchmarks/rtg_vcfeval/{caller}/`):
- `tp.vcf.gz` — True Positives (callset)
- `tp-baseline.vcf.gz` — True Positives (baseline)
- `fn.vcf.gz` — False Negatives
- `fp.vcf.gz` — False Positives
- `summary.txt` — Precision, Recall, F-measure

#### 5.2 Hợp nhất các file output fn/fp/tp-callset/tp-baseline sau khi dùng rtgtool/vcfeval của 4 công cụ gọi biến thể

``` bash
# pwd: variant-benchmarking

OUTDIR=results/benchmarks/rtg_vcfeval/merged
mkdir -p "$OUTDIR"
```

```bash
# pwd: variant-benchmarking
# FN

FN_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*fn*.vcf.gz | head -n 1)
FN_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*fn*.vcf.gz | head -n 1)
FN_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*fn*.vcf.gz | head -n 1)
FN_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*fn*.vcf.gz | head -n 1)

for X in gatk:"$FN_GATK" deepvariant:"$FN_DV" strelka2:"$FN_ST" freebayes:"$FN_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.fn.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.fn.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/fn_4callers.vcf.gz" \
  "$OUTDIR/gatk.fn.vcf.gz" "$OUTDIR/deepvariant.fn.vcf.gz" \
  "$OUTDIR/strelka2.fn.vcf.gz" "$OUTDIR/freebayes.fn.vcf.gz"
tabix -f -p vcf "$OUTDIR/fn_4callers.vcf.gz"
```

```bash
# pwd: variant-benchmarking
# FP

FP_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*fp*.vcf.gz | head -n 1)
FP_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*fp*.vcf.gz | head -n 1)
FP_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*fp*.vcf.gz | head -n 1)
FP_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*fp*.vcf.gz | head -n 1)

for X in gatk:"$FP_GATK" deepvariant:"$FP_DV" strelka2:"$FP_ST" freebayes:"$FP_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.fp.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.fp.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/fp_4callers.vcf.gz" \
  "$OUTDIR/gatk.fp.vcf.gz" "$OUTDIR/deepvariant.fp.vcf.gz" \
  "$OUTDIR/strelka2.fp.vcf.gz" "$OUTDIR/freebayes.fp.vcf.gz"
tabix -f -p vcf "$OUTDIR/fp_4callers.vcf.gz"
```

```bash
# pwd: variant-benchmarking
# TP_CALLSET (TP này là TP_callset, không phải baseline)

TP_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*tp*.vcf.gz | grep -v baseline | head -n 1)
TP_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*tp*.vcf.gz | grep -v baseline | head -n 1)
TP_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*tp*.vcf.gz | grep -v baseline | head -n 1)
TP_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*tp*.vcf.gz | grep -v baseline | head -n 1)

for X in gatk:"$TP_GATK" deepvariant:"$TP_DV" strelka2:"$TP_ST" freebayes:"$TP_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.tp_callset.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.tp_callset.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/tp_callset_4callers.vcf.gz" \
  "$OUTDIR/gatk.tp_callset.vcf.gz" "$OUTDIR/deepvariant.tp_callset.vcf.gz" \
  "$OUTDIR/strelka2.tp_callset.vcf.gz" "$OUTDIR/freebayes.tp_callset.vcf.gz"
tabix -f -p vcf "$OUTDIR/tp_callset_4callers.vcf.gz"
```

```bash
# pwd: variant-benchmarking
# TP_BASELINE

TPB_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*tp-baseline*.vcf.gz | head -n 1)
TPB_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*tp-baseline*.vcf.gz | head -n 1)
TPB_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*tp-baseline*.vcf.gz | head -n 1)
TPB_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*tp-baseline*.vcf.gz | head -n 1)

for X in gatk:"$TPB_GATK" deepvariant:"$TPB_DV" strelka2:"$TPB_ST" freebayes:"$TPB_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.tp_baseline.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.tp_baseline.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/tp_baseline_4callers.vcf.gz" \
  "$OUTDIR/gatk.tp_baseline.vcf.gz" "$OUTDIR/deepvariant.tp_baseline.vcf.gz" \
  "$OUTDIR/strelka2.tp_baseline.vcf.gz" "$OUTDIR/freebayes.tp_baseline.vcf.gz"
tabix -f -p vcf "$OUTDIR/tp_baseline_4callers.vcf.gz"
```

Output (`results/benchmarks/rtg_vcfeval/merged/`):
- `fn_4callers.vcf.gz` — FN hợp nhất 4 caller
- `fp_4callers.vcf.gz` — FP hợp nhất 4 caller
- `tp_callset_4callers.vcf.gz` — TP callset hợp nhất
- `tp_baseline_4callers.vcf.gz` — TP baseline hợp nhất

### 5.3 Chuyển file VCF sang CSV để phục vụ cho công đoạn xử lý dữ liệu và trực quan hoá

```bash
# pwd: variant-benchmarking
# FN

IN=results/benchmarks/rtg_vcfeval/merged/fn_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/fn_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```

```bash
# pwd: variant-benchmarking
# FP

IN=results/benchmarks/rtg_vcfeval/merged/fp_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/fp_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```

```bash
# pwd: variant-benchmarking
# TP_CALLSET

IN=results/benchmarks/rtg_vcfeval/merged/tp_callset_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/tp_callset_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```


```bash
# pwd: variant-benchmarking
# TP_BASELINE

IN=results/benchmarks/rtg_vcfeval/merged/tp_baseline_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/tp_baseline_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```

Output (`results/benchmarks/rtg_vcfeval/merged/`):
- `fn_4callers.csv` — FN dạng CSV
- `fp_4callers.csv` — FP dạng CSV
- `tp_callset_4callers.csv` — TP callset dạng CSV
- `tp_baseline_4callers.csv` — TP baseline dạng CSV

### Phần VI. Đánh giá mức độ nguy hiểm của lỗi FN/FP (Functional Risk of Errors)

Bên cạnh các metric cổ điển (Precision/Recall/F1), phần này bổ sung đánh giá
**"ai sai nguy hiểm hơn"** thông qua 2 lớp phân tích:

- **Lớp 1 — SnpEff + ACMG**: phân loại lỗi theo impact (HIGH/MODERATE/LOW/MODIFIER) và phân loại lâm sàng (P/LP/VUS/LB/B)
- **Lớp 2 — AlphaGenome**: chấm điểm tác động chức năng bằng deep learning, đặc biệt mạnh với biến thể novel và non-coding

#### 6.1. Chuẩn bị cơ sở dữ liệu cho SnpEff / SnpSift

```bash
# pwd: variant-benchmarking

# Download database GRCh38
java -jar snpEff.jar download GRCh38.mane.1.0.ensembl

# Download dbNSFP (SIFT, PolyPhen2, PROVEAN, MutationTaster, CADD, MetaLR)
wget -P data/snpeff/data https://dbnsfp.s3.amazonaws.com/dbNSFP4.9a.zip

# Trích chr22 từ file zip (header + dữ liệu chr22), không cần giải nén toàn bộ
unzip -p data/snpeff/data/dbNSFP4.9a.zip dbNSFP4.9a_variant.chr22.gz \
  | zcat | head -1 | bgzip > data/snpeff/data/dbNSFP4.9a.chr22.txt.gz
unzip -p data/snpeff/data/dbNSFP4.9a.zip dbNSFP4.9a_variant.chr22.gz \
  | zcat | sed 1d | bgzip >> data/snpeff/data/dbNSFP4.9a.chr22.txt.gz
tabix -s 1 -b 2 -e 2 data/snpeff/data/dbNSFP4.9a.chr22.txt.gz

# Xoá file zip gốc (~35GB)
rm data/snpeff/data/dbNSFP4.9a.zip

# Download ClinVar
mkdir -p data/snpeff/data/clinvar
wget -P data/snpeff/data/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -P data/snpeff/data/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

#### 6.2. Chú giải biến thể bằng SnpEff + SnpSift

Chạy pipeline chú giải cho **toàn bộ** FN/FP của từng caller:
SnpEff (chức năng) → dbSNP (rsID) → dbNSFP (SIFT, PolyPhen-2, PROVEAN, MutationTaster, CADD, MetaLR, 1000Gp3_AF) → ClinVar.

```bash
# pwd: variant-benchmarking
bash stages/07_snpeff_annotate.sh
```

Output: `results/annotation/{caller}/fn.annotated.vcf.gz` và `fp.annotated.vcf.gz`

#### 6.3. Phân tích Functional Risk (Lớp 1: SnpEff Impact + ACMG)

Phân tích toàn bộ FN/FP (không lọc trước) theo 2 chiều:
- **SnpEff Impact Distribution**: đếm biến thể HIGH / MODERATE / LOW / MODIFIER per caller
- **ACMG 5-tier Classification**: phân loại P / LP / VUS / LB / B per caller

Đánh giá chức năng bằng 6 công cụ: SIFT (≤ 0.05), PolyPhen-2 (≥ 0.85), PROVEAN (≤ -2.5), MutationTaster (D/A), CADD (≥ 20), MetaLR (≥ 0.5).

```bash
# pwd: variant-benchmarking
python3 scripts/functional_risk_analysis.py
```

Output:
- `results/annotation/impact_distribution.csv` — phân bố impact per caller
- `results/annotation/acmg_summary.csv` — thống kê ACMG per caller
- `results/annotation/{caller}/{fn|fp}.acmg.tsv` — chi tiết từng biến thể

#### 6.4. Chấm điểm AlphaGenome (Lớp 2: Deep Learning Functional Impact)

AlphaGenome bổ sung cho SnpEff ở những nơi SnpEff yếu: biến thể **novel** và **non-coding/regulatory**. Chạy riêng cho FN/FP của từng caller.

```bash
# pwd: variant-benchmarking
# Chạy cho từng caller × category (fn/fp)
python3 scripts/batch_alphagenome_fn_fp.py --caller gatk --error-type fn
python3 scripts/batch_alphagenome_fn_fp.py --caller gatk --error-type fp
python3 scripts/batch_alphagenome_fn_fp.py --caller deepvariant --error-type fn
python3 scripts/batch_alphagenome_fn_fp.py --caller deepvariant --error-type fp
python3 scripts/batch_alphagenome_fn_fp.py --caller strelka2 --error-type fn
python3 scripts/batch_alphagenome_fn_fp.py --caller strelka2 --error-type fp
python3 scripts/batch_alphagenome_fn_fp.py --caller freebayes --error-type fn
python3 scripts/batch_alphagenome_fn_fp.py --caller freebayes --error-type fp
```

Output: `results/alphagenome/{caller}_{fn|fp}_scores.csv`

