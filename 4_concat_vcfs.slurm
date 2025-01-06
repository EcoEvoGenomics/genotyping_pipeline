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

# User variables
PIPELINE_REPOSITORY_DIR=/path/to/pipeline/repository
VCF_LIST=/path/to/list/of/vcf/files/to/concatenate
OUTPUT_VCF=vcf_concat/NAME_HERE.vcf.gz
OUTPUT_VCF_NORM=vcf_concat/NAME_HERE_norm.vcf.gz

# Prepare environment
set -o errexit
set -o nounset
module --quiet purge

# Load modules
module load BCFtools/1.19-GCC-13.2.0
module list

# Work in Nextflow pipeline directory
cd $PIPELINE_REPOSITORY_DIR

if [ ! -e vcf_concat/ ]
then
    mkdir vcf_concat
fi

# Concatenate and index VCFs
bcftools concat -f $VCF_LIST --threads 8 -n -O z -o $OUTPUT_VCF
bcftools index $OUTPUT_VCF

# Next, normalise
bcftools norm -d none -O z -o $OUTPUT_VCF_NORM $OUTPUT_VCF
bcftools index $OUTPUT_VCF_NORM
