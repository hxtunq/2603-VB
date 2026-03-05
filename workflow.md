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
mkdir -p results/variants
mkdir -p results/plots
mkdir -p logs
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
rm -f dbsnp138.hg38.vcf.gz dbsnp138.hg38.vcf.gz.tbi
rm -f Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
rm -f 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

```bash
# pwd: variant-benchmarking/data/reference

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

```bash
# pwd: variant-benchmarking/data/reference

# Tạo high-confidence BED = non-N regions trừ đi vùng low-mappability
# Dùng làm confident region (-f) cho hap.py
bedtools subtract \
    -a chr22_non_N_regions.bed \
    -b umap_k100_mappability.chr22.bed.gz \
    > chr22_highconf.bed
```

Output:
- `data/reference/dbsnp138.hg38.chr22.vcf.gz` (+`.tbi`) — dbSNP chr22
- `data/reference/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz` (+`.tbi`) — Mills indels chr22
- `data/reference/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz` (+`.tbi`) — 1000G SNPs chr22
- `data/reference/chr22_non_N_regions.bed` — vùng non-N
- `data/reference/chr22_highconf.bed` — high-confidence region (non-N − low-mappability) cho hap.py

#### 2.3. Tạo đột biến SNPs và INDELs với simuG

Thông số được chỉnh theo dữ liệu thực tế (chr22 ~50 Mb, ~20K variants cho GIAB samples).

```bash
# pwd: variant-benchmarking/data

# clone simuG về folder <data>
git clone https://github.com/yjx1217/simuG.git

# dùng simuG tạo đột biến — thông số realistic
perl simuG/simuG.pl \
  -refseq reference/chr22.fa \
  -snp_count 18000 \
  -titv_ratio 2.0 \
  -indel_count 2500 \
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

#### 2.5. Download BED files cho exome targets và stratification

Lấy từ AstraZeneca-NGS reference_data (hg38), extract chr22.
**Phải chạy trước ART-Illumina** vì WES simulation cần exome BED để extract vùng target.

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

#### 2.6. Tạo file fastq giả lập bằng ART Illumina — WGS

Tạo reads WGS từ **toàn bộ** mutated genome (chr22 ~50Mb):
- **WGS coverage:** 10x, 20x, 30x, 50x (tương tự GIAB WGS 22–37x)

```bash
# pwd: variant-benchmarking/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

# WGS: coverage thấp-trung bình, chạy trên toàn bộ chr22
for COV in 10 20 30 50; do
  echo "=== WGS: Generating ${COV}x coverage ==="

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

#### 2.7. Tạo file fastq giả lập bằng ART Illumina — WES

Tạo reads WES chỉ từ **vùng exome** (~1-2Mb trên chr22):
- **WES coverage:** 100x, 200x (tương tự GIAB WES 183–249x)

> **Lưu ý:** Dùng `bedtools getfasta` để extract chỉ vùng exome từ mutated genome
> trước khi cho vào ART. Nhờ đó 200x chỉ cover ~1-2Mb thay vì ~50Mb → tiết kiệm
> ~25-50x dung lượng và thời gian. Reads sau khi align (BWA-MEM) vẫn map lên toàn
> bộ chr22, nhưng coverage chỉ cao ở vùng exome → đúng pattern WES thực tế.

```bash
# pwd: variant-benchmarking/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"
EXOME_BED="reference/Exome-Agilent_V6.chr22.bed"

# Index mutated genome cho bedtools
samtools faidx "${MUTATED_FASTA}"

# Extract chỉ vùng exome từ mutated genome
bedtools getfasta \
    -fi "${MUTATED_FASTA}" \
    -bed "${EXOME_BED}" \
    -fo "${SIM_DIR}/${PREFIX}_exome.fa"

# WES: coverage cao, chỉ trên vùng exome
for COV in 50 100 200; do
  echo "=== WES: Generating ${COV}x coverage ==="

  art_illumina \
      -ss HS25 \
      -i "${SIM_DIR}/${PREFIX}_exome.fa" \
      -p \
      -l 150 \
      -f ${COV} \
      -m 350 \
      -s 50 \
      -rs 42 \
      -o "${SIM_DIR}/${PREFIX}_${COV}x_wes_" \
      -na

  # đổi tên và nén
  mv "${SIM_DIR}/${PREFIX}_${COV}x_wes_1.fq" "${SIM_DIR}/${PREFIX}_${COV}x_wes_R1.fastq"
  mv "${SIM_DIR}/${PREFIX}_${COV}x_wes_2.fq" "${SIM_DIR}/${PREFIX}_${COV}x_wes_R2.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_wes_R1.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_wes_R2.fastq"

  echo "Done: ${SIM_DIR}/${PREFIX}_${COV}x_wes_R{1,2}.fastq.gz"
done
```

Output:
- `data/simulated/SIMULATED_SAMPLE_chr22_exome.fa` — mutated genome chỉ vùng exome
- `data/simulated/SIMULATED_SAMPLE_chr22_{COV}x_wes_R1.fastq.gz` — WES paired-end read 1
- `data/simulated/SIMULATED_SAMPLE_chr22_{COV}x_wes_R2.fastq.gz` — WES paired-end read 2

### Phần III. Tiền xử lý dữ liệu cho công đoạn gọi biến thể — Multi-coverage

Theo paper (Barbitoff et al. 2022), preprocessing chung cho **tất cả callers** chỉ gồm:
**BWA-MEM → sort → MarkDuplicates → coverage stats**

> **Lưu ý quan trọng:**  
> **BQSR** (Base Quality Score Recalibration) chỉ áp dụng cho **GATK HaplotypeCaller** (theo GATK Best Practices).  
> Các caller khác (DeepVariant, Strelka2, FreeBayes, DNAscope) dùng trực tiếp **dedup BAM** mà không cần BQSR.  
> BQSR đã được tích hợp sẵn trong script `pipelines/03_call_hc.sh`.

Chạy lặp cho **mỗi mức coverage**, chia thành WGS và WES:

#### 3.1 BWA-MEM alignment + sort

```bash
# pwd: variant-benchmarking

PREFIX="SIMULATED_SAMPLE_chr22"
REF="data/reference/chr22.fa"

# === WGS (10x, 20x, 30x, 50x) ===
for COV in 10 20 30 50; do
  echo "=== BWA-MEM WGS: ${COV}x ==="

  R1="data/simulated/${PREFIX}_${COV}x_R1.fastq.gz"
  R2="data/simulated/${PREFIX}_${COV}x_R2.fastq.gz"
  OUTDIR="results/preprocessing/${COV}x"
  mkdir -p "${OUTDIR}" "logs/${COV}x"

  bwa mem -t 4 -M \
    -R "@RG\tID:${PREFIX}_${COV}x\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${R1}" "${R2}" \
    2> "logs/${COV}x/bwa_mem.log" | \
    samtools sort -@ 4 -m 2G -o "${OUTDIR}/${PREFIX}_aligned.bam" -

  samtools index "${OUTDIR}/${PREFIX}_aligned.bam"
done

# === WES (100x, 200x) — reads từ vùng exome, align lên toàn bộ chr22 ===
for COV in 50 100 200; do
  echo "=== BWA-MEM WES: ${COV}x ==="

  R1="data/simulated/${PREFIX}_${COV}x_wes_R1.fastq.gz"
  R2="data/simulated/${PREFIX}_${COV}x_wes_R2.fastq.gz"
  OUTDIR="results/preprocessing/${COV}x_wes"
  mkdir -p "${OUTDIR}" "logs/${COV}x_wes"

  bwa mem -t 4 -M \
    -R "@RG\tID:${PREFIX}_${COV}x_wes\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${R1}" "${R2}" \
    2> "logs/${COV}x_wes/bwa_mem.log" | \
    samtools sort -@ 4 -m 2G -o "${OUTDIR}/${PREFIX}_aligned.bam" -

  samtools index "${OUTDIR}/${PREFIX}_aligned.bam"
done
```

#### 3.2 MarkDuplicates

```bash
# pwd: variant-benchmarking

PREFIX="SIMULATED_SAMPLE_chr22"

# === WGS ===
for COV in 10 20 30 50; do
  echo "=== MarkDuplicates WGS: ${COV}x ==="
  OUTDIR="results/preprocessing/${COV}x"

  gatk MarkDuplicates \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_dedup.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true \
    2>&1 | tee "logs/${COV}x/markduplicates.log"
done

# === WES ===
for COV in 50 100 200; do
  echo "=== MarkDuplicates WES: ${COV}x ==="
  OUTDIR="results/preprocessing/${COV}x_wes"

  gatk MarkDuplicates \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_dedup.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true \
    2>&1 | tee "logs/${COV}x_wes/markduplicates.log"
done
```

#### 3.3 Coverage statistics

```bash
# pwd: variant-benchmarking

PREFIX="SIMULATED_SAMPLE_chr22"

# === WGS ===
for COV in 10 20 30 50; do
  echo "=== Stats WGS: ${COV}x ==="
  OUTDIR="results/preprocessing/${COV}x"

  samtools stats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
  samtools flagstat "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
  samtools idxstats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_idxstats.txt"

  mosdepth -t 4 --by 1000 \
    "${OUTDIR}/${PREFIX}_coverage" \
    "${OUTDIR}/${PREFIX}_dedup.bam"

  echo "FINAL_BAM=${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/bam_path.sh"
done

# === WES ===
for COV in 50 100 200; do
  echo "=== Stats WES: ${COV}x ==="
  OUTDIR="results/preprocessing/${COV}x_wes"

  samtools stats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
  samtools flagstat "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
  samtools idxstats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_idxstats.txt"

  mosdepth -t 4 --by 1000 \
    "${OUTDIR}/${PREFIX}_coverage" \
    "${OUTDIR}/${PREFIX}_dedup.bam"

  echo "FINAL_BAM=${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/bam_path.sh"
done
```

> Preprocessing đã được mô tả chi tiết trong các bước 3.1–3.3 ở trên.

Output:
- WGS: `results/preprocessing/{COV}x/` — BAM, stats, coverage (COV = 10, 20, 30, 50)
- WES: `results/preprocessing/{COV}x_wes/` — BAM, stats, coverage (COV = 100, 200)

### Phần IV. Gọi biến thể — Multi-coverage, WGS + WES

Chạy tất cả caller cho mỗi mức coverage:
- **WGS:** 10x, 20x, 30x, 50x — dùng BAM từ `results/preprocessing/{COV}x/`
- **WES:** 100x, 200x — dùng BAM từ `results/preprocessing/{COV}x_wes/`

#### 4.1 WGS Variant Calling

```bash
# pwd: variant-benchmarking

for COV in 10 20 30 50; do
  echo "=== Variant Calling: ${COV}x WGS ==="

  bash pipelines/03_call_hc.sh ${COV} > "logs/${COV}x/gatk.log" 2>&1              # GATK HC (3 filters)
  bash pipelines/04_call_dv.sh ${COV} > "logs/${COV}x/deepvariant.log" 2>&1               # DeepVariant
  bash pipelines/05_call_strelka.sh ${COV} > "logs/${COV}x/strelka2.log" 2>&1           # Strelka2 (with Manta)
  bash pipelines/06_call_freebayes.sh ${COV} > "logs/${COV}x/freebayes.log" 2>&1         # FreeBayes
  bash pipelines/07_call_dnascope.sh ${COV} > "logs/${COV}x/dnascope.log" 2>&1          # DNAscope (optional)
done
```

#### 4.2 WES Variant Calling (100x, 200x)

```bash
# pwd: variant-benchmarking

for COV in 50 100 200; do
  echo "=== Variant Calling: ${COV}x WES ==="

  bash pipelines/03_call_hc_wes.sh ${COV} > "logs/${COV}x_wes/gatk_wes.log" 2>&1          # GATK HC WES
  bash pipelines/04_call_dv_wes.sh ${COV} > "logs/${COV}x_wes/deepvariant_wes.log" 2>&1           # DeepVariant WES
  bash pipelines/05_call_strelka_wes.sh ${COV} > "logs/${COV}x_wes/strelka2_wes.log" 2>&1      # Strelka2 WES (with Manta --exome)
  bash pipelines/06_call_freebayes_wes.sh ${COV} > "logs/${COV}x_wes/freebayes_wes.log" 2>&1    # FreeBayes WES
  bash pipelines/07_call_dnascope_wes.sh ${COV} > "logs/${COV}x_wes/dnascope_wes.log" 2>&1     # DNAscope WES (optional)
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
> **Total VCFs:** 4 coverages × WGS (5) + 2 coverages × WES (5) = 30+ VCFs

### Phần V. Đánh giá kết quả gọi biến thể bằng hap.py

Sử dụng [hap.py](https://github.com/Illumina/hap.py) (via Docker) để so sánh kết quả variant calling với truth VCF.
hap.py tích hợp sẵn vcfeval engine và tính Precision, Recall, F1 cho cả SNP và INDEL.

- **WGS:** đánh giá trên toàn bộ high-confidence region (`-f chr22_highconf.bed`)
- **WES:** thu hẹp thêm vào exome targets (`-T Exome-Agilent_V6.chr22.bed`)
- **Stratification:** phân tầng kết quả theo CDS, Exome, Low-mappability (nếu có `stratification_chr22.tsv`)

#### 5.1 Chạy hap.py cho tất cả callers × coverages

```bash
# pwd: variant-benchmarking
bash evaluation/eval_happy.sh
```

Script tự động lặp qua:
- WGS: 10x, 20x, 30x, 50x × 5 callers (gatk, deepvariant, strelka2, freebayes, dnascope)
- WES: 100x, 200x × 5 callers (gatk_wes, deepvariant_wes, strelka2_wes, freebayes_wes, dnascope_wes)

Output per caller: `results/eval/{COV}x/{WGS|WES}/{caller}_eval_data/report.*`

#### 5.2 Tổng hợp thống kê vào file TSV

```bash
# pwd: variant-benchmarking
bash evaluation/gather_stats.sh
```

Output: `results/eval/all_stats.tsv` — bảng tổng hợp với các cột: Coverage, Mode, Caller, + toàn bộ metrics từ hap.py report

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

