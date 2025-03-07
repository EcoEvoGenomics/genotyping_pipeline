#!/bin/bash
# written 29/06/2023

# nextflow management script 

# Job name:
#SBATCH --job-name=india
#
# Project:
#SBATCH --account=nn10082k
#
# Wall clock limit (hh:mm:ss):
#SBATCH --time=120:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1

 # notify end of job
#SBATCH --mail-user=mark.ravinet@ibv.uio.no
#SBATCH --mail-type=FAIL

module load Miniconda3/22.11.1-1
export CONDA_PKGS_DIRS=/cluster/projects/nn10082k/conda_users/msravine/package-cache
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/nextflow

cd /path/to/working/directory/

REF=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa

# run nextflow pipeline - not recommended to run all in one go - comment out commands you don't want to run
nextflow run 1_trim_map_realign.nf --samples /path/to/samples.csv --ref $REF 
#nextflow run 2_call_variants.nf --bams genotyping_cram_list.txt --ref $REF --windows sparrow_genome_windows.list
#nextflow run 3_filter_variants.nf  --miss 0.5 --min_depth 3 --max_depth 30
