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

    # Path of repository directory
    repository_path=

    # Path to .CSV with each row as "ID, F_READ_PATH, R_READ_PATH"
    sample_csv=

    # Which steps to run? Note: combine_vcfs is optional, but steps must be run in order.
    trim_align=yes
    call_vcf=yes
    filt_vcf=yes
    combine_vcfs=yes
    multiqc=yes

    # Reference genome
    ref_genome=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
    ref_index=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa.fai
    ref_scaffold_name='scaffold'
    ref_ploidy_file=./pipeline/defaults/default.ploidy

    # Threshold depth of CRAMs?
    downsample_large_crams=no
    max_cram_depth=30

    # Genotyping settings
    window_size=10000000

    # Filtering settings
    vcf_filt_miss=0.8
    vcf_filt_q_site_ps=30
    vcf_filt_q_site_gs=30
    vcf_filt_min_depth=5
    vcf_filt_max_depth=30
    vcf_filt_min_geno_depth=5
    vcf_filt_max_geno_depth=30
    vcf_filt_keep=""
    stats_downsample_sites=10000
    
    # NB: You should not change these unless you know what you are doing.
    output_dir=./output
    trim_align_output_dir=${output_dir}/trimmed_and_aligned
    call_vcf_output_dir=${output_dir}/called_variants

    # You can change this, though, if you want to refilter your VCF!
    filt_vcf_output_dir=${output_dir}/filtered_variants

### --- End user input --- ###

# Prepare environment
set -o errexit
set -o nounset
module --quiet purge

# Load modules
module load Miniconda3/22.11.1-1

# Activate conda environment
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_users/eriksro/genotyping_pipeline_experimental

# Function to handle missing output directories
mkmissingdir() {
    if [ ! -e $1 ]; then
        mkdir $1
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

    printf "miss %s\nq_site_ps %s\nq_site_gs %s\nmin_depth %s\nmax_depth %s\nmin_geno_depth %s\nmax_geno_depth %s\nkeep %s\n" \
        "$vcf_filt_miss" \
        "$vcf_filt_q_site_ps" \
        "$vcf_filt_q_site_gs" \
        "$vcf_filt_min_depth" \
        "$vcf_filt_max_depth" \
        "$vcf_filt_min_geno_depth" \
        "$vcf_filt_max_geno_depth" \
        "$vcf_filt_keep" \
        > $filt_vcf_output_dir/filt_params.txt

    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/filter_variants.nf \
        -c ./pipeline/config/filter_variants.config \
        --vcf_dir $call_vcf_output_dir \
        --miss $vcf_filt_miss \
        --q_site_ps $vcf_filt_q_site_ps \
        --q_site_gs $vcf_filt_q_site_gs \
        --min_depth $vcf_filt_min_depth \
        --max_depth $vcf_filt_max_depth \
        --min_geno_depth $vcf_filt_min_geno_depth \
        --max_geno_depth $vcf_filt_max_geno_depth \
        --keep $vcf_filt_keep \
        --publish_dir $filt_vcf_output_dir \
        --stats_downsample_sites $stats_downsample_sites
fi

if [ $combine_vcfs = 'yes' ]; then
    chkprevious "Step: combine_vcfs" $filt_vcf_output_dir

    nextflow -log ./.nextflow/nextflow.log \
        run ./pipeline/nextflow/combine_vcf.nf \
        -c ./pipeline/config/combine_vcf.config \
        --input_dir $filt_vcf_output_dir/vcf_filtered \
        --pop_structure_label '_ps' \
        --genome_scan_label '_gs' \
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
