# Variant Calling Benchmark

Hệ thống benchmark so sánh hiệu suất các công cụ gọi biến thể (variant caller) trên dữ liệu giả lập từ chromosome 22 (hg38). Thiết kế theo phương pháp luận của [Barbitoff et al. 2022, BMC Genomics](https://doi.org/10.1186/s12864-022-08365-3).

## Tổng quan

```
chr22.fa (hg38)
    │
    ├── simuG → 18K SNPs + 2.5K INDELs (truth VCF)
    │
    └── ART Illumina (6 mức coverage: 10x, 20x, 30x, 50x, 100x, 200x)
            │
            └── BWA-MEM → MarkDuplicates → dedup BAM
                    │
                    ├── GATK HC: BQSR → HaplotypeCaller → CNN/Hard Filter
                    ├── DeepVariant, Strelka2+Manta, FreeBayes, DNAscope
                    │
                    ├── WGS: 5 callers × 6 coverages
                    └── WES: 5 callers × 3 coverages (≥50x)
                            │
                            └── hap.py + vcfeval → F1, Precision, Recall
```

**Tổng cộng:** ~45 VCFs → đánh giá cả WGS lẫn WES với stratification regions.

## Cấu trúc thư mục

```
variant-calling-benchmark/
├── config/
│   └── config.sh                    # Cấu hình chung (resources, paths, params)
│
├── pipelines/                       # Scripts gọi biến thể
│   ├── 03_call_hc.sh                # GATK HC (BQSR + Hard Filter)
│   ├── 03_call_hc_wes.sh            # GATK HC — WES mode (-L exome BED)
│   ├── 04_call_dv.sh                # DeepVariant WGS
│   ├── 04_call_dv_wes.sh            # DeepVariant WES (--model_type=WES)
│   ├── 05_call_strelka.sh           # Manta + Strelka2 WGS
│   ├── 05_call_strelka_wes.sh       # Manta + Strelka2 WES (--exome)
│   ├── 06_call_freebayes.sh         # FreeBayes WGS
│   ├── 06_call_freebayes_wes.sh     # FreeBayes WES (--targets)
│   ├── 07_call_dnascope.sh          # Sentieon DNAscope WGS (optional)
│   └── 07_call_dnascope_wes.sh      # Sentieon DNAscope WES (optional)
│
├── evaluation/
│   ├── eval_happy.sh                # hap.py evaluation (multi-coverage, WGS+WES)
│   └── gather_stats.sh              # Thu thập kết quả vào TSV
│
├── visualization/
│   ├── benchmark_plots.R            # R script vẽ biểu đồ benchmark
│   └── plot_summary.py              # Python summary plots
│
├── workflow.md                      # Quy trình chi tiết từng bước (copy-paste ready)
└── caller_benchmark-main/           # Reference code từ Barbitoff et al. 2022
```

## Variant Callers

| # | Caller | WGS | WES | Filter strategy |
|---|--------|-----|-----|-----------------|
| 1 | **GATK HaplotypeCaller** | ✅ | ✅ | Hard Filter (GATK Best Practices) |
| 2 | **DeepVariant** v1.9.0 | ✅ | ✅ | Default (model-based) |
| 3 | **Strelka2** v2.9.10 | ✅ | ✅ | Default (with Manta indel candidates) |
| 4 | **FreeBayes** v1.3.10 | ✅ | ✅ | QUAL + strand bias filters |
| 5 | **DNAscope** (Sentieon) | ✅ | ✅ | ML model (optional, cần license) |

## Coverage Design

| Coverage | WGS eval | WES eval | Mục đích |
|----------|----------|----------|----------|
| 10x | ✅ | — | Low coverage stress test |
| 20x | ✅ | — | Gần GIAB WGS (~22x) |
| 30x | ✅ | — | Standard clinical WGS |
| 50x | ✅ | ✅ | High WGS / Low WES |
| 100x | ✅ | ✅ | Medium WES |
| 200x | ✅ | ✅ | Gần GIAB WES (~200x) |

## Cách chạy

Xem chi tiết trong [workflow.md](workflow.md). Tóm tắt:

```bash
# 1. Setup
source config/config.sh

# 2. Simulate (simuG → ART multi-coverage)
# Xem workflow.md §2.3–2.5

# 3. Preprocessing (per coverage)
# Xem workflow.md Part III

# 4. Variant calling
for COV in 10 20 30 50 100 200; do
  export PREPROC_DIR="results/preprocessing/${COV}x"
  source "${PREPROC_DIR}/bam_path.sh"
  bash pipelines/03_call_hc.sh
  bash pipelines/04_call_dv.sh
  bash pipelines/05_call_strelka.sh
  bash pipelines/06_call_freebayes.sh
  bash pipelines/07_call_dnascope.sh  # optional
done

# 5. Evaluate
bash evaluation/eval_happy.sh

# 6. Visualize
Rscript visualization/benchmark_plots.R
```

## Đánh giá

- **Tool:** hap.py + RTGtools vcfeval
- **Metrics:** F1 Score, Precision, Recall (SNP + INDEL riêng)
- **Stratification:** CDS regions, Exome targets, Low-mappability regions
- **Statistical power:** Nhiều coverage levels → boxplot + so sánh xu hướng

## Tham khảo

```bash
# WGS coverages
for cov in 10 20 30 50; do
    bash pipelines/03_call_hc.sh $cov
done

# WES coverages
for cov in 50 100 200; do
    bash pipelines/03_call_hc_wes.sh $cov
done
```
