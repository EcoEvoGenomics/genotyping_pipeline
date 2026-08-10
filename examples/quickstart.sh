#!/bin/bash

# ADMIN
#SBATCH --job-name=XENO
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k

# RESOURCE ALLOCATION
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=99:00:00

# This quickstart script works on the NRIS Saga HPC for users with access
# to the nn10082k project number. Must be run from top-level in the repo.

module --quiet purge
module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/Nextflow25.04.6
bash ./XENO
