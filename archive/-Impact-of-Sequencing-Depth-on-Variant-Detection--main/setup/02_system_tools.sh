#!/bin/bash
set -euo pipefail

# Make conda available in non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "Installing system tools..."
sudo apt-get update
sudo apt-get install -y \
  wget \
  git \
  google-cloud-sdk

# NOTE: tools requiring ngs1 environment
conda activate ngs1

echo "Downloading Picard (required for uBAM filtering)..."
wget -q https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar \
  -O ~/picard.jar

echo "Installing FreeBayes..."
conda install -y -c bioconda freebayes
conda deactivate

# NOTE: tools requiring cromwell environment (Java 11)
conda activate cromwell

echo "Downloading Cromwell jar..."
wget -q https://github.com/broadinstitute/cromwell/releases/download/81/cromwell-81.jar \
  -O cromwell.jar
conda deactivate

echo "Pulling DeepVariant Docker image..."
docker pull google/deepvariant:latest

echo "Tool setup completed successfully"

