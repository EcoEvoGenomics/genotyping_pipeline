#!/bin/bash

# Admin:
#SBATCH --job-name=concat_vcf
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k

# Resource allocation:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12g
#SBATCH --time=02:00:00

# Prepare environment
set -o errexit
set -o nounset
module --quiet purge

# Load modules
module load BCFtools/1.19-GCC-13.2.0
module list

# Work in Nextflow pipeline directory
cd $PIPELINE_REPOSITORY_DIR

# Concatenate and index VCFs
bcftools concat -f $VCF_LIST --threads 8 -n -O z -o $OUTPUT_VCF
bcftools index $OUTPUT_VCF

# Next, normalise
bcftools norm -d none -O z -o $OUTPUT_VCF_NORM $OUTPUT_VCF
bcftools index $OUTPUT_VCF_NORM
