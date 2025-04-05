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
#SBATCH --time=99:00:00

### ----- User input ----- ###

    # Path to .CSV of input samples
    # Note: samples are rows and the columns "ID, F_READ_PATH, R_READ_PATH" without headers
    sample_csv=

    # Which steps to run?
    # Note: steps must be run in order, but you can repeat filt_vcf and multiqc with different settings
    trim_align=yes
    call_vcf=yes
    filt_vcf=yes
    multiqc=yes

    # Settings for step trim_align
    # Note: max_cram_depth only applies if downsample_large_crams=yes
    downsample_large_crams=no
    max_cram_depth=30

    # Settings for step call_vcf
    window_size=10000000

    # Settings for step filt_vcf
    # Note: change filtering_label and re-run filt_vcf to refilter output from call_vcf
    filtering_label='DEFAULT_POP_STRUCTURE'
    filtering_min_alleles=2
    filtering_max_alleles=2
    filtering_miss=0.8
    filtering_q_site=30
    filtering_min_depth=5
    filtering_max_depth=30
    filtering_min_geno_depth=5
    filtering_max_geno_depth=30
    filtering_keep=''

    # Reference genome
    ref_genome=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
    ref_index=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa.fai
    ref_scaffold_name='scaffold'
    ref_ploidy_file=./pipeline/defaults/default.ploidy
    
    # Paths to this repository and to a suitable conda environment
    repository_path=
    conda_environment_path=/cluster/projects/nn10082k/conda_users/eriksro/genotyping_pipeline_experimental

### --- End user input --- ###

# Prepare environment
set -o errexit
set -o nounset
module --quiet purge

# Load modules
module load Miniconda3/22.11.1-1

# Activate conda environment
source ${EBROOTMINICONDA3}/bin/activate
conda activate ${conda_environment_path}

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

output_dir=./output
trim_align_output_dir=${output_dir}/01-aligned_reads
call_vcf_output_dir=${output_dir}/02-variants_unfiltered
filt_vcf_output_dir=${output_dir}/03-variants_filtered/${filtering_label}

mkmissingdir $output_dir

if [ $trim_align = 'yes' ]; then
    mkmissingdir $trim_align_output_dir

    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/trim_and_align.nf \
        -c ./pipeline/config/trim_and_align.config \
        --ref $ref_genome \
        --ref_scaffold_name $ref_scaffold_name \
        --samples $sample_csv \
        --downsample_crams $downsample_large_crams \
        --max_cram_depth $max_cram_depth \
        --publish_dir $trim_align_output_dir
fi

if [ $call_vcf = 'yes' ]; then
    chkprevious "Step: call_vcf" $trim_align_output_dir
    mkmissingdir $call_vcf_output_dir

    windows_dir=$call_vcf_output_dir/genome_windows
    mkmissingdir $windows_dir
    bash ./pipeline/shell/create_genome_windows.sh $ref_index $window_size $ref_scaffold_name $windows_dir

    cram_list=${call_vcf_output_dir}/genotyped_crams.list
    find $PWD/$trim_align_output_dir -name '*.cram' > $cram_list

    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/call_variants.nf \
        -c ./pipeline/config/call_variants.config \
        --ref $ref_genome \
        --crams $cram_list \
        --windows_dir $windows_dir \
        --ploidy_file $ref_ploidy_file \
        --publish_dir $call_vcf_output_dir
fi

if [ $filt_vcf = 'yes' ]; then
    chkprevious "Step: filt_vcf" $call_vcf_output_dir
    mkmissingdir $filt_vcf_output_dir

    printf "miss %s\nq_site %s\nmin_depth %s\nmax_depth %s\nmin_geno_depth %s\nmax_geno_depth %s\nkeep %s\n" \
        "$filtering_miss" \
        "$filtering_q_site" \
        "$filtering_min_depth" \
        "$filtering_max_depth" \
        "$filtering_min_geno_depth" \
        "$filtering_max_geno_depth" \
        "$filtering_keep" \
        > $filt_vcf_output_dir/filt_params.txt

    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/filter_variants.nf \
        -c ./pipeline/config/filter_variants.config \
        --vcf_dir $call_vcf_output_dir/chroms \
        --filtering_label $filtering_label \
        --min_alleles= $filtering_min_alleles \
        --max_alleles= $filtering_max_alleles \
        --miss $filtering_miss \
        --q_site $filtering_q_site \
        --min_depth $filtering_min_depth \
        --max_depth $filtering_max_depth \
        --min_geno_depth $filtering_min_geno_depth \
        --max_geno_depth $filtering_max_geno_depth \
        --keep $filtering_keep \
        --publish_dir $filt_vcf_output_dir
fi

if [ $multiqc = 'yes' ]; then
    # NB: Temporary solution only!
    # Deactivate pipeline environment, then conda. Purge and use Saga MultiQC module to get later version.
    conda deactivate
    conda deactivate
    module --quiet purge
    module load MultiQC/1.22.3-foss-2023b
    multiqc --outdir $output_dir $output_dir
fi

# End work
