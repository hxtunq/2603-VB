# Variant Calling Performance Analysis: Impact of Sequencing Depth on SNP and Indel Detection

## Overview

How much sequencing depth is really enough for reliable variant calling?

This project provides a comprehensive benchmarking study using the GIAB HG002 exome (Chr22) to evaluate how sequencing depth impacts SNP and indel detection accuracy, while comparing three widely used variant callers: **GATK HaplotypeCaller**, **DeepVariant**, and **FreeBayes**.

A key technical aspect of this work was adapting GATK Best Practices pipelines originally designed for cloud environments to run efficiently on a local 16 GB RAM laptop. This included:

* Converting legacy GRCh37-aligned BAMs to uBAM and realigning to GRCh38 to eliminate reference and alignment bias
* Optimizing Cromwell/WDL workflows, Docker memory limits, and Java heap sizes for local execution under strict resource constraints
* Maintaining full reproducibility and benchmarking accuracy despite non-cloud hardware

## Key Findings

* **Ultra-low coverage (≤10×)** leads to severe sensitivity loss and unreliable genotyping
* **~40× coverage** represents an optimal balance for most research applications
* **≥80× coverage** is required for robust indel detection and clinical-grade confidence
* **DeepVariant** consistently achieved the highest precision, with dramatically fewer false positives than FreeBayes and slightly outperforming GATK

## Methodology

The project was benchmarked against GIAB truth sets using hap.py, with orthogonal quality assessment (Ts/Tv, FP/FN rates), and executed end-to-end using Docker, Conda, GATK Best Practices, and reproducible workflows.

## Repository Contents

* **`Final_report_ghada_ahmed.pdf`** - Complete analysis report with detailed methodology, results, visualizations, and discussion
* **`setup/`** - Configuration files and environment setup for reproducible execution
* **`scripts/`** - Analysis pipelines, variant calling workflows, and benchmarking code
* **`csv_summaries/`** - Performance metrics, summary statistics, and tabular results

## Technologies Used

- GATK Best Practices (Cromwell/WDL)
- DeepVariant, GATK HaplotypeCaller, FreeBayes
- hap.py benchmarking toolkit
- Docker, Conda
- GRCh38 reference genome
- GIAB HG002 v4.2.1 truth set

## Quick Start

For detailed setup instructions, analysis workflows, and complete results, please refer to the final report and contents of the `setup/` and `scripts/` directories.

---

**Author:** Ghada Ahmed  
**Project Type:** Bioinformatics Benchmarking Study  
**Date:** jan-2026
