#!/bin/bash

# Create conda environment
mkdir envs
env=envs/default
conda create --prefix $env -y

# Activate environment and install packages
source activate $env

conda install -c conda-forge r-base=4.1.2 -y
conda install -c conda-forge r-seurat=4.1.0 -y
conda install -c bioconda bioconductor-batchelor=1.10.0 -y

# Create hidden file to indicate active environment
echo $env > .env
