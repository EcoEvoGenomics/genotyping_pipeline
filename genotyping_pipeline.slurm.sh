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
#SBATCH --mem-per-cpu=4G
#SBATCH --time=99:00:00

### SETTINGS (1 / 2)
###
### ----------------- User input ----------------- ###

    # Path to .CSV of input samples
    # Note: samples are rows and the columns "ID, LANE, F_READ_PATH, R_READ_PATH" without headers
    sample_csv=
    
    # Path to this repository
    repository_path=

    # Which steps to run?
    # Note: steps must be run in order, but you can repeat filter_variants and multiqc with different settings
    trim_align_reads=yes
    call_variants=yes
    filter_variants=yes
    multiqc=yes

    # Settings for step trim_align_reads
    # Note: read_target only applies if downsample_reads=yes
    # Note: flag 0x400 is for optical and PCR duplicates
    deduplicate_reads=no
    downsample_reads=no
    read_target=1000000
    exclude_flags=0x400

    # Settings for step call_variants
    window_size=10000000
    concatenate_unfiltered_vcfs=no

    # Settings for step filter_variants
    # Note: change filtering_label and re-run filter_variants to refilter output from call_variants
    filtering_label='default_filters'
    filtering_min_alleles=2
    filtering_max_alleles=2
    filtering_max_missing=0.8
    filtering_min_meanDP=5
    filtering_max_meanDP=30
    filtering_minDP=5
    filtering_maxDP=30
    filtering_minQ=30
    filtering_keep=''

    # Reference genome
    ref_genome=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
    ref_index=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa.fai
    ref_scaffold_name='scaffold'
    ref_ploidy_file=./pipeline/assets/default.ploidy
    
    # Path to a conda environment with Nextflow
    nextflow_envir_path=/cluster/projects/nn10082k/conda_group/nextflow

    # Name of Nextflow configuration profile to use (see nextflow.config)
    nextflow_profile="saga"

### --------------- End user input --------------- ###

### SETTINGS (2 / 2)
###
### ------------ Set up environment -------------- ###

    # This script must be run from an environment with:
    # - Conda
    # - Singularity (default on Saga)
    # - Slurm (default on Saga)

    module --quiet purge
    module load Miniconda3/22.11.1-1
    source ${EBROOTMINICONDA3}/bin/activate

### ------------ End set up environment ---------- ###

# Prepare environment
set -o errexit
set -o nounset
conda activate ${nextflow_envir_path}

# Function to handle missing output directories
mkmissingdir() {
    if [ ! -e $1 ]; then
        mkdir -p $1
    fi
}

# Function to check previous step is done
chkprevious() {
    if [ ! -e $2 ]; then
        echo "Error. ${1} expected output from previous step to exist in directory ${2}."
        exit
    fi
}

# Begin work
cd $repository_path

output_dir=${repository_path}/output
trim_align_output_dir=${output_dir}/01-aligned_reads
call_variants_output_dir=${output_dir}/02-variants_unfiltered
filter_variants_output_dir=${output_dir}/03-variants_filtered/${filtering_label}
multiqc_output_dir=${output_dir}/04-multiqc

mkmissingdir $output_dir

if [ $trim_align_reads = 'yes' ]; then
    mkmissingdir $trim_align_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/trim_align_reads.nf \
        -with-report $trim_align_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        --samples $sample_csv \
        --deduplicate $deduplicate_reads \
        --downsample $downsample_reads \
        --read_target $read_target \
        --exclude_flags $exclude_flags \
        --ref_genome $ref_genome \
        --ref_index $ref_index \
        --ref_scaffold_name $ref_scaffold_name \
        --publish_dir $trim_align_output_dir
fi

if [ $call_variants = 'yes' ]; then
    chkprevious "Step: call_variants" $trim_align_output_dir
    mkmissingdir $call_variants_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/call_variants.nf \
        -with-report $call_variants_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        --cram_dir $trim_align_output_dir \
        --window_size $window_size \
        --ref_genome $ref_genome \
        --ref_index $ref_index \
        --ref_scaffold_name $ref_scaffold_name \
        --ref_ploidy_file $ref_ploidy_file \
        --concatenate_vcf $concatenate_unfiltered_vcfs \
        --publish_dir $call_variants_output_dir
fi

if [ $filter_variants = 'yes' ]; then
    chkprevious "Step: filter_variants" $call_variants_output_dir
    mkmissingdir $filter_variants_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/filter_variants.nf \
        -with-report $filter_variants_output_dir/workflow_report.html \
        -profile $nextflow_profile \
        --vcf_dir $call_variants_output_dir/chroms \
        --filtering_label $filtering_label \
        --min_alleles $filtering_min_alleles \
        --max_alleles $filtering_max_alleles \
        --max_missing $filtering_max_missing \
        --min_meanDP $filtering_min_meanDP \
        --max_meanDP $filtering_max_meanDP \
        --minDP $filtering_minDP \
        --maxDP $filtering_maxDP \
        --minQ $filtering_minQ \
        --keep $filtering_keep \
        --publish_dir $filter_variants_output_dir
fi

if [ $multiqc = 'yes' ]; then
    mkmissingdir $multiqc_output_dir
    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/run_multiqc.nf \
        -profile $nextflow_profile \
        --results_dir $output_dir \
        --publish_dir $multiqc_output_dir
fi

# End work
