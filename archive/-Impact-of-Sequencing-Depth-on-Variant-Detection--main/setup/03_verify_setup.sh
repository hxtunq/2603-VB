#!/bin/bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

echo "Checking system setup..."

echo "Checking conda"
conda --version

echo "Checking ngs1 environment"
conda activate ngs1
samtools --version
freebayes --version
conda deactivate

echo "Checking hapy_py27 environment"
conda activate hapy_py27
hap.py --help | head -n 1
conda deactivate

echo "Checking cromwell environment"
conda activate cromwell
java -version
conda deactivate

echo "All checks passed"
