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

    # Path to .CSV with each row as "ID, F_READ_PATH, R_READ_PATH, ADAPTER"
    sample_csv=

    # Which steps to run?
    trim_align=yes
    downsample_bams=yes
    call_vcf=yes
    filt_vcf=yes
    concat_vcf=no # This step currently does not work

    # Reference genome
    ref_genome=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
    ref_index=/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa.fai
    ref_scaffold_name='scaffold'
    ref_ploidy_file=./pipeline/defaults/default.ploidy

    # Directory of trimming adapters
    adapter_dir=/cluster/projects/nn10082k/trimmomatic_adapters/

    # BAM maximum depth before downsampling
    bam_max_depth=20

    # Genotyping settings
    window_size=10000000

    # Filtering settings
    vcf_filt_miss=0.8
    vcf_filt_q_site1=30
    vcf_filt_q_site2=30
    vcf_filt_min_depth=5
    vcf_filt_max_depth=30
    vcf_filt_min_geno_depth=5
    vcf_filt_max_geno_depth=30
    vcf_filt_keep=""
    
    # NB: Only change if you know what you are doing (ask Erik), typically for re-filtering a VCF
    #     Also, . must be the genotyping_pipeline repository
    output_dir=./output
    trim_align_output_dir=${output_dir}/trim_align
    call_vcf_output_dir=${output_dir}/raw_vcf
    filt_vcf_output_dir=${output_dir}/filt_vcf
    concat_vcf_output_dir=${filt_vcf_output_dir}

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

# Begin work
cd $repository_path
mkmissingdir $output_dir

if [ $trim_align = 'yes' ]; then
    mkmissingdir $trim_align_output_dir
    nextflow run ./pipeline/nextflow/trim_and_align.nf \
        -c ./pipeline/config/trim_and_align.config \
        --ref $ref_genome \
        --samples $sample_csv \
        --trim $adapter_dir \
        --publish_dir $trim_align_output_dir
fi

if [ $downsample_bams = 'yes' ]; then
    
    if [ ! -e $trim_align_output_dir ]; then
        echo "Error. Expected to downsample BAMs, but no BAMs exist in specified directory."
        exit 1
    fi

    undownsampled_bams=${trim_align_output_dir}/undownsampled_bams.list
    find $PWD/$trim_align_output_dir/align/ -name '*.*am' > $undownsampled_bams

    nextflow run ./pipeline/nextflow/downsample_bams.nf \
        -c ./pipeline/config/downsample_bams.config \
        --depth $bam_max_depth \
        --bams $undownsampled_bams \
        --publish_dir $trim_align_output_dir
fi

if [ $call_vcf = 'yes' ]; then
    mkmissingdir $call_vcf_output_dir
    windows_dir=$call_vcf_output_dir/genome_windows
    mkmissingdir $windows_dir
    bash ./pipeline/shell/create_genome_windows.sh $ref_index $window_size $ref_scaffold_name $windows_dir

    bam_list=${call_vcf_output_dir}/genotyped_bams.list

    if [ -e ${trim_align_output_dir}/downsample_align ]; then
        echo "Downsampled BAMs exist in ${trim_align_output_dir}/downsample_align. Genotyping downsampled BAMs ..."
        find $PWD/$trim_align_output_dir/align_downsample/ -name '*.*am' > $bam_list
    else
        echo "Genotyping non-downsampled BAMs from ${trim_align_output_dir} ..."
        find $PWD/$trim_align_output_dir/align/ -name '*.*am' > $bam_list
    fi

    nextflow run ./pipeline/nextflow/call_variants.nf \
        -c ./pipeline/config/call_variants.config \
        --ref $ref_genome \
        --bams $bam_list \
        --windows_dir $windows_dir \
        --ploidyFile $ref_ploidy_file \
        --publish_dir $call_vcf_output_dir
fi

if [ $filt_vcf = 'yes' ]; then
    mkmissingdir $filt_vcf_output_dir
    echo "VCF filtering parameters:\n" \
        "miss" $vcf_filt_miss "\n" \
        "q_site1" $vcf_filt_q_site1 "\n" \
        "q_site2" $vcf_filt_q_site2 "\n" \
        "min_depth" $vcf_filt_min_depth "\n" \
        "max_depth" $vcf_filt_max_depth "\n" \
        "min_geno_depth" $vcf_filt_min_geno_depth "\n" \
        "max_geno_depth" $vcf_filt_max_geno_depth "\n" \
        "keep" $vcf_filt_keep "\n" \
        > $filt_vcf_output_dir/filt_params.txt
    nextflow run ./pipeline/nextflow/filter_variants.nf \
        -c ./pipeline/config/filter_variants.config \
        --vcf_dir $call_vcf_output_dir/vcf \
        --miss $vcf_filt_miss \
        --q_site1 $vcf_filt_q_site1 \
        --q_site2 $vcf_filt_q_site2 \
        --min_depth $vcf_filt_min_depth \
        --max_depth $vcf_filt_max_depth \
        --min_geno_depth $vcf_filt_min_geno_depth \
        --max_geno_depth $vcf_filt_max_geno_depth \
        --keep $vcf_filt_keep \
        --publish_dir $filt_vcf_output_dir
fi

if [ $concat_vcf = 'yes' ]; then
    mkmissingdir $concat_vcf_output_dir
    vcf_list=${concat_vcf_output_dir}/concatenated_vcfs.list
    find $PWD/$filt_vcf_output_dir/vcf_filtered/ -name '*.vcf' > $vcf_list
    sbatch ./pipeline/shell/concat_vcfs.slurm.sh \
        --export=\
        PIPELINE_REPOSITORY_DIR=$repository_path,\
        VCF_LIST=$vcf_list,\
        OUTPUT_VCF=$concat_vcf_output_dir/concatenated.vcf.gz,\
        OUTPUT_VCF_NORM=$concat_vcf_output_dir/concatenated_normalised.vcf.gz
fi

# End work
