#!/bin/bash

#SBATCH --account=its
#SBATCH --cpus-per-task=24
#SBATCH --gpus=1

spack load snakemake \
singularityce \
miniconda3 \
graphviz \
trf \
samtools@1.16.1 \
cdhit \
gffread \
bedops \
minimap2^python@3.10 \
py-biopython \
py-torch/rvl \
py-pybedtools \
py-scikit-learn \
py-statsmodels

snakemake --use-conda --conda-frontend conda --cores "$SLURM_CPUS_PER_TASK" "$@"

# Run example:
# sbatch --time=07-00:00:00 run_centromer_detection.sh \
#   results/Guy11_chr1/CENTROMERE_SCORING/Guy11_chr1_1000/centro_candidates.bed


