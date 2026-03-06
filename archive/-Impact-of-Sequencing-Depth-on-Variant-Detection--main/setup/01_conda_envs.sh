#!/bin/bash
set -euo pipefail

# Initialize conda for non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "=== Creating conda environments (run once per machine) ==="

# Core NGS tools
conda create -y -n ngs1 \
  python=3.10 \
  samtools \
  bcftools \
  bedtools \
  htslib

# GIAB benchmarking (hap.py)
conda create -y -n hapy_py27 \
  python=2.7 \
  -c bioconda \
  hap.py

# Cromwell / WDL execution
conda create -y -n cromwell \
  -c conda-forge \
  openjdk=11

echo "=== Conda environments created successfully ==="
