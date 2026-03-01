# Variant Calling Benchmark

So sánh và đánh giá hiệu năng của 4 công cụ gọi biến thể — **GATK 4.6.2.0**, **DeepVariant 1.8.0**, **Strelka2 2.9.10**, **FreeBayes 1.3.10** — trên dữ liệu WGS giả lập (chr22, hg38, 60x).

## Cấu trúc thư mục

```text
variant-calling-benchmark/
│
├── readme.md                          # Tổng quan project
├── workflow.md                        # Quy trình thực hiện chi tiết (step-by-step)
│
├── config/
│   └── config.sh                      # Cấu hình chung: paths, versions, tham số simulation & calling
│
├── stages/                            # Các script chạy variant calling
│   ├── 03_variant_calling_gatk.sh     # GATK HaplotypeCaller + hard filtering
│   ├── 04_variant_calling_deepvariant.sh  # DeepVariant via Docker
│   ├── 05_variant_calling_strelka2.sh # Strelka2 Germline via Docker
│   ├── 06_variant_calling_freebayes.sh    # FreeBayes
│   └── 07_snpeff_annotate.sh          # SnpEff + SnpSift annotation (FN/FP)
│
├── scripts/                           # Scripts phụ trợ
│   ├── helper_functions.sh            # Logging, timer, resource tracking
│   ├── normalize_vcf.sh               # Chuẩn hoá VCF cho benchmarking (bcftools norm)
│   ├── rename_chromosomes.sh          # Đổi tên chromosome (22 → chr22) nếu cần
│   ├── functional_risk_analysis.py    # Phân tích risk: SnpEff Impact + ACMG 5-tier
│   └── batch_alphagenome_fn_fp.py     # Chấm điểm FN/FP bằng AlphaGenome (Deep Learning)
│
├── data-viz/                          # Trực quan hoá dữ liệu
│   ├── upset-plot.Rmd                 # UpSet plot (R Markdown)
│   ├── upset_long_baseline.csv        # Dữ liệu TP-baseline cho UpSet plot
│   └── upset_long_callset.csv         # Dữ liệu TP-callset cho UpSet plot
│
├── logs/                              # Log từ các bước chạy
│   └── runtime.csv                    # Thời gian chạy của mỗi caller (giây)
│
├── data/                              # ⬇️ Sinh ra khi chạy workflow
│   ├── reference/                     # Reference genome (chr22.fa), index, known-sites (dbSNP, Mills, 1000G)
│   └── simulated/                     # Truth VCF + FASTQ giả lập từ simuG & ART
│
└── results/                           # ⬇️ Kết quả sau khi chạy pipeline
    ├── preprocessing/                 # BAM sau alignment + BQSR, QC metrics
    │   ├── fastqc_raw/                # FastQC report cho raw reads
    │   ├── *_dup_metrics.txt          # Thống kê duplicate reads
    │   ├── *_recal.table              # Bảng BQSR recalibration
    │   ├── *_stats.txt                # samtools stats
    │   ├── *_flagstat.txt             # samtools flagstat
    │   ├── *_coverage.*               # Coverage analysis (mosdepth)
    │   └── bam_path.sh                # Export đường dẫn BAM cuối cho các caller
    │
    ├── variants/                      # VCF output của từng caller
    │   ├── gatk/                      # *_gatk_pass.norm.vcf.gz
    │   ├── deepvariant/               # *_deepvariant_pass.norm.vcf.gz
    │   ├── strelka2/                  # *_strelka2_pass.norm.vcf.gz
    │   └── freebayes/                 # *_freebayes_pass.norm.vcf.gz
    │
    └── benchmarks/
        └── rtg_vcfeval/               # So sánh từng caller vs truth (RTG VCFeval)
            ├── gatk/                  # TP/FP/FN + summary.txt
            ├── deepvariant/           # TP/FP/FN + summary.txt
            ├── strelka2/              # TP/FP/FN + summary.txt
            ├── freebayes/             # TP/FP/FN + summary.txt
            └── merged/                # FN/FP/TP hợp nhất 4 caller (VCF + CSV)
```

## Quy trình

Toàn bộ quy trình được mô tả chi tiết trong [`workflow.md`](workflow.md), gồm 6 phần:

1. **Tạo cấu trúc folder**
2. **Download & chuẩn bị dữ liệu** — reference genome, known sites, simuG, ART
3. **Tiền xử lý** — FastQC → BWA-MEM → MarkDuplicates → BQSR → Filter MQ≥20
4. **Gọi biến thể** — GATK, DeepVariant, Strelka2, FreeBayes
5. **Benchmarking** — RTG VCFeval (TP/FP/FN, Precision/Recall/F1)
6. **Đánh giá Functional Risk** — SnpEff + ACMG + AlphaGenome
