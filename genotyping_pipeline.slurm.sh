#!/bin/bash

# ADMIN
#SBATCH --job-name=genotyping
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k
#
# RESOURCE ALLOCATION
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00

# User definitions
pipeline_directory=
sample_csv=
output_dir=

# Prepare environment
set -o errexit
set -o nounset
module --quiet purge

# Load modules
module load Miniconda3/22.11.1-1

# Activate conda environment
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/nextflow

# Begin work
cd $pipeline_directory

# Trim and align
nextflow run ./pipeline/nextflow/trim_and_align.nf \
    -c ./pipeline/config/trim_and_align.config \
    --samples $sample_csv \
    --publish_dir $output_dir

# Call variants
nextflow run ./pipeline/nextflow/call_variants.nf \
    -c ./pipeline/config/call_variants.config \
    --bams $bams \
    --windows $windows

# Filter variants
nextflow run ./pipeline/nextflow/filter_variants.nf \
    -c ./pipeline/config/filter_variants.nf \
    --filtering-params $filtering_params

# End work
