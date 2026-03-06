## Workflow

This file documents the implemented benchmark only. Annotation-heavy stages such as SnpEff, dbNSFP, ClinVar, and AlphaGenome are future work and are not part of the runnable pipeline in this repo.

## 1. Setup

```bash
git clone <repo-url> variant-calling-benchmark
cd variant-calling-benchmark

mkdir -p data/reference
mkdir -p data/simulated
mkdir -p results/preprocessing
mkdir -p results/variants
mkdir -p results/eval
mkdir -p results/eval_track_b
mkdir -p results/plots
mkdir -p logs

source config/config.sh
```

## 2. Reference Assets

### 2.1 Reference FASTA

```bash
cd data/reference

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz

samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

### 2.2 Known Sites For BQSR

```bash
cd data/reference

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.vcf.gz
bcftools view -r chr22 -Oz -o dbsnp138.hg38.chr22.vcf.gz dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.chr22.vcf.gz

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
bcftools view -r chr22 -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
bcftools view -r chr22 -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz \
  1000G_phase1.snps.high_confidence.hg38.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
```

### 2.3 BED Assets And hap.py Support Files

```bash
cd data/reference

wget -O CDS-canonical.bed \
  https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/CDS-canonical.bed
grep -w "chr22" CDS-canonical.bed > CDS-canonical.chr22.bed
rm CDS-canonical.bed

wget -O Exome-Agilent_V6.bed \
  https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/Exome-Agilent_V6.bed
grep -w "chr22" Exome-Agilent_V6.bed > Exome-Agilent_V6.chr22.bed
rm Exome-Agilent_V6.bed

wget -O umap_k100_mappability.bed.gz \
  https://github.com/AstraZeneca-NGS/reference_data/raw/master/hg38/tricky_regions/umap_k100_mappability.bed.gz
wget -O umap_k100_mappability.bed.gz.tbi \
  https://github.com/AstraZeneca-NGS/reference_data/raw/master/hg38/tricky_regions/umap_k100_mappability.bed.gz.tbi
tabix umap_k100_mappability.bed.gz chr22 | bgzip > umap_k100_mappability.chr22.bed.gz
tabix -p bed umap_k100_mappability.chr22.bed.gz

python3 - <<'PY'
from Bio import SeqIO
import re

with open("chr22_highconf_tmp.bed", "w") as out:
    for record in SeqIO.parse("chr22.fa", "fasta"):
        for match in re.finditer(r"[ACGT]+", str(record.seq).upper()):
            if match.end() - match.start() >= 1000:
                out.write(f"{record.id}\t{match.start()}\t{match.end()}\n")
PY

bedtools subtract -a chr22_highconf_tmp.bed -b umap_k100_mappability.chr22.bed.gz > chr22_highconf.bed
rm chr22_highconf_tmp.bed

cat > stratification_chr22.tsv <<'EOF'
CDS	CDS-canonical.chr22.bed
Exome_Targets	Exome-Agilent_V6.chr22.bed
Low_Mappability	umap_k100_mappability.chr22.bed.gz
EOF

rtg format -o chr22.sdf chr22.fa
```

`chr22.sdf` is required before running hap.py with `--engine vcfeval`.

### 2.4 Sentieon Assets

Track A DNAscope and Track B pangenome both require a valid license and model bundles.

```bash
export SENTIEON_LICENSE=your_server:port

# Download and unpack the official model bundles into data/reference/models/
# - DNAscopeIlluminaWGS2.0.bundle
# - DNAscopeIlluminaWES2.0.bundle

# Build the pangenome index prefix used by the supplementary Track B scripts.
# The population VCF must match the reference build and contig naming.
sentieon bwa-mem2-pangenome index data/reference/chr22.fa population_variants.vcf.gz
```

The resulting pangenome prefix should match `PANGENOME_INDEX` in `config/config.sh`.

## 3. Truth Set And Simulated FASTQs

### 3.1 Truth VCF From simuG

```bash
cd data

git clone https://github.com/yjx1217/simuG.git

perl simuG/simuG.pl \
  -refseq reference/chr22.fa \
  -snp_count 18000 \
  -indel_count 2500 \
  -prefix simulated/SIMULATED_SAMPLE_chr22

samtools faidx reference/chr22.fa
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.vcf
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.vcf

bcftools concat \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf | \
bcftools sort -Oz -o simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

tabix -p vcf simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz
```

### 3.2 WGS FASTQs

```bash
cd data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

for COV in 10 20 30 50; do
  art_illumina \
    -ss HS25 \
    -i "${MUTATED_FASTA}" \
    -p \
    -l 150 \
    -f "${COV}" \
    -m 350 \
    -s 50 \
    -rs 42 \
    -o "${SIM_DIR}/${PREFIX}_${COV}x_" \
    -na

  mv "${SIM_DIR}/${PREFIX}_${COV}x_1.fq" "${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq"
  mv "${SIM_DIR}/${PREFIX}_${COV}x_2.fq" "${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_R1.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_R2.fastq"
done
```

### 3.3 WES FASTQs

```bash
cd data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"
EXOME_BED="reference/Exome-Agilent_V6.chr22.bed"

samtools faidx "${MUTATED_FASTA}"
bedtools getfasta -fi "${MUTATED_FASTA}" -bed "${EXOME_BED}" -fo "${SIM_DIR}/${PREFIX}_exome.fa"

for COV in 50 100 200; do
  art_illumina \
    -ss HS25 \
    -i "${SIM_DIR}/${PREFIX}_exome.fa" \
    -p \
    -l 150 \
    -f "${COV}" \
    -m 350 \
    -s 50 \
    -rs 42 \
    -o "${SIM_DIR}/${PREFIX}_${COV}x_wes_" \
    -na

  mv "${SIM_DIR}/${PREFIX}_${COV}x_wes_1.fq" "${SIM_DIR}/${PREFIX}_${COV}x_wes_R1.fastq"
  mv "${SIM_DIR}/${PREFIX}_${COV}x_wes_2.fq" "${SIM_DIR}/${PREFIX}_${COV}x_wes_R2.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_wes_R1.fastq"
  gzip "${SIM_DIR}/${PREFIX}_${COV}x_wes_R2.fastq"
done
```

## 4. Shared Preprocessing

The shared preprocessing BAM is the only valid input for Track A comparisons.

### 4.1 Alignment And Sort

```bash
PREFIX="SIMULATED_SAMPLE_chr22"
REF="data/reference/chr22.fa"

for COV in 10 20 30 50; do
  R1="data/simulated/${PREFIX}_${COV}x_R1.fastq.gz"
  R2="data/simulated/${PREFIX}_${COV}x_R2.fastq.gz"
  OUTDIR="results/preprocessing/${COV}x"
  mkdir -p "${OUTDIR}" "logs/${COV}x"

  bwa mem -t 4 -M \
    -R "@RG\tID:${PREFIX}_${COV}x\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${R1}" "${R2}" | \
    samtools sort -@ 4 -m 2G -o "${OUTDIR}/${PREFIX}_aligned.bam" -

  samtools index "${OUTDIR}/${PREFIX}_aligned.bam"
done

for COV in 50 100 200; do
  R1="data/simulated/${PREFIX}_${COV}x_wes_R1.fastq.gz"
  R2="data/simulated/${PREFIX}_${COV}x_wes_R2.fastq.gz"
  OUTDIR="results/preprocessing/${COV}x_wes"
  mkdir -p "${OUTDIR}" "logs/${COV}x_wes"

  bwa mem -t 4 -M \
    -R "@RG\tID:${PREFIX}_${COV}x_wes\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${REF}" "${R1}" "${R2}" | \
    samtools sort -@ 4 -m 2G -o "${OUTDIR}/${PREFIX}_aligned.bam" -

  samtools index "${OUTDIR}/${PREFIX}_aligned.bam"
done
```

### 4.2 MarkDuplicates And Coverage Stats

```bash
PREFIX="SIMULATED_SAMPLE_chr22"

for COV in 10 20 30 50; do
  OUTDIR="results/preprocessing/${COV}x"

  gatk MarkDuplicates \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_dedup.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true

  samtools stats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
  samtools flagstat "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
  samtools idxstats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_idxstats.txt"
  mosdepth -t 4 --by 1000 "${OUTDIR}/${PREFIX}_coverage" "${OUTDIR}/${PREFIX}_dedup.bam"

  echo "FINAL_BAM=${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/bam_path.sh"
done

for COV in 50 100 200; do
  OUTDIR="results/preprocessing/${COV}x_wes"

  gatk MarkDuplicates \
    --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
    -I "${OUTDIR}/${PREFIX}_aligned.bam" \
    -O "${OUTDIR}/${PREFIX}_dedup.bam" \
    -M "${OUTDIR}/${PREFIX}_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true

  samtools stats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_stats.txt"
  samtools flagstat "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_flagstat.txt"
  samtools idxstats "${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/${PREFIX}_idxstats.txt"
  mosdepth -t 4 --by 1000 "${OUTDIR}/${PREFIX}_coverage" "${OUTDIR}/${PREFIX}_dedup.bam"

  echo "FINAL_BAM=${OUTDIR}/${PREFIX}_dedup.bam" > "${OUTDIR}/bam_path.sh"
done
```

## 5. Variant Calling

### 5.1 Track A: Shared BAM

```bash
for COV in 10 20 30 50; do
  bash pipelines/03_call_hc.sh "${COV}"
  bash pipelines/04_call_dv.sh "${COV}"
  bash pipelines/05_call_strelka.sh "${COV}"
  bash pipelines/06_call_freebayes.sh "${COV}"
  bash pipelines/07_call_dnascope.sh "${COV}"
done

for COV in 50 100 200; do
  bash pipelines/03_call_hc_wes.sh "${COV}"
  bash pipelines/04_call_dv_wes.sh "${COV}"
  bash pipelines/05_call_strelka_wes.sh "${COV}"
  bash pipelines/06_call_freebayes_wes.sh "${COV}"
  bash pipelines/07_call_dnascope_wes.sh "${COV}"
done
```

### 5.2 Track B: Supplementary DNAscope Pangenome

```bash
for COV in 10 20 30 50; do
  bash pipelines/07_call_dnascope_pangenome.sh "${COV}"
done

for COV in 50 100 200; do
  bash pipelines/07_call_dnascope_pangenome_wes.sh "${COV}"
done
```

Track B must not be pooled into the main caller ranking because it uses its own alignment pipeline from FASTQ.

## 6. Evaluation

### 6.1 Main Track

```bash
bash evaluation/eval_happy.sh
```

Outputs:

- Evaluation tree: `results/eval/...`
- Aggregated stats: `results/eval/all_stats.tsv`

### 6.2 Supplementary Track

```bash
bash evaluation/eval_happy_track_b.sh
```

Outputs:

- Evaluation tree: `results/eval_track_b/...`
- Aggregated stats: `results/eval_track_b/all_stats.tsv`

`evaluation/gather_stats.sh` is the single collector and can also be called manually:

```bash
bash evaluation/gather_stats.sh results/eval results/eval/all_stats.tsv
bash evaluation/gather_stats.sh results/eval_track_b results/eval_track_b/all_stats.tsv
```

## 7. Visualization

```bash
Rscript visualization/benchmark_plots.R results/eval/all_stats.tsv results/plots
python visualization/plot_summary.py results/eval/all_stats.tsv results/plots
```

Both plotting scripts now create separate output sets per `Mode + Coverage` combination to avoid mixing WGS and WES or collapsing multiple coverages into the same figure.

## 8. Method Notes

- Track A is the paper-grade benchmark because all callers receive the same dedup BAM.
- Track B is useful to show Sentieon’s end-to-end pangenome best case, but it is not a fair caller-only comparison.
- The archived HG002 depth study in `archive/-Impact-of-Sequencing-Depth-on-Variant-Detection--main/` is exome-only chr22 downsampling. It supports qualitative depth intuition only:
  - `<=10x` is unstable
  - `~40x` is workable for general research use
  - indels benefit from higher depth

## 9. Future Work

Planned but not implemented in this repo:

- SnpEff / SnpSift annotation
- dbNSFP / ClinVar integration
- ACMG-style FN/FP risk summaries
- AlphaGenome scoring
