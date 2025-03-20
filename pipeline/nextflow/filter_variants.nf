#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Filter VCF
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default parameters
params.vcf_dir = "${params.publish_dir}/call_vcf/vcf"
params.publish_dir = './output'

// Default filtering parameters
params.miss=0.8
params.q_site1=30
params.q_site2=30
params.min_depth=5
params.max_depth=30
params.min_geno_depth=5
params.max_geno_depth=30
params.keep=""

// Workflow
workflow{
  Channel.fromPath("${params.vcf_dir}/*.vcf.gz") | rm_indels | vcf_filter
}

// Filtering, Step 1 - Remove spanning indels and re-normalise
process rm_indels {

  input:
  path (anno_vcf)

  output:
  tuple \
    file ("${anno_vcf.simpleName}.vcf.gz"), \
    file ("${anno_vcf.simpleName}.vcf.gz.csi") 

  script:
  """
  bcftools view --threads ${task.cpus} -V indels -e 'ALT="*" | N_ALT>1' $anno_vcf | bcftools norm --threads ${task.cpus} -D -O z -o ${anno_vcf.simpleName}.vcf.gz

  bcftools index --threads ${task.cpus} ${anno_vcf.simpleName}.vcf.gz
  """
}

// Filtering, Step 2 - Filtering for population structure analyses and for genome scans
process vcf_filter {

  publishDir "${params.publish_dir}/vcf_filtered", saveAs: { filename -> "$filename" }, mode: 'copy'

  input:
  tuple file(rm_indel_vcf), file(rm_indel_vcf_index) 

  output:
  tuple \
    file ("${rm_indel_vcf.simpleName}_filtered_ps.vcf.gz"), \
    file ("${rm_indel_vcf.simpleName}_filtered_ps.vcf.gz.csi"), \
    file ("${rm_indel_vcf.simpleName}_filtered_gs.vcf.gz"), \
    file ("${rm_indel_vcf.simpleName}_filtered_gs.vcf.gz.csi") 

  script:
  """
  if [[ -f ${params.keep} ]]; then

    echo "File of individuals to filter provided - adding --keep option."

    # for pop structure
    vcftools --gzvcf $rm_indel_vcf  --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing ${params.miss} \
    --minQ ${params.q_site1} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filtered_ps.vcf.gz

    bcftools index --threads ${task.cpus} ${rm_indel_vcf.simpleName}_filtered_ps.vcf.gz

    # for genome scans
    vcftools --gzvcf $rm_indel_vcf --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --max-alleles 2 \
    --minQ ${params.q_site2} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filtered_gs.vcf.gz
    
    bcftools index --threads ${task.cpus} ${rm_indel_vcf.simpleName}_filtered_gs.vcf.gz

  else

    echo "Not filtering for specific individuals"
    # for pop structure
    vcftools --gzvcf $rm_indel_vcf  --remove-indels --remove-filtered-all \
    --min-alleles 2 --max-alleles 2 \
    --max-missing ${params.miss} --minQ ${params.q_site1} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filtered_ps.vcf.gz

    bcftools index --threads ${task.cpus} ${rm_indel_vcf.simpleName}_filtered_ps.vcf.gz

    # for genome scans
      vcftools --gzvcf $rm_indel_vcf --remove-indels --remove-filtered-all \
    --max-alleles 2 \
    --minQ ${params.q_site2} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filtered_gs.vcf.gz
    
    bcftools index --threads ${task.cpus} ${rm_indel_vcf.simpleName}_filtered_gs.vcf.gz

  fi
  """
}

// Run QC for unfiltered VCFs
process qc_pre {

  input:
  tuple file(vcf), file(csi)
  path windows_dir
  each window

  output:
  path "${vcf.simpleName}_${window}_qc_pre.txt"

  script:
  """
  bcftools stats --threads ${task.cpus} --regions ${window} ${vcf} > ${vcf.simpleName}_${window}_qc_pre.txt
  """
}

// Run QC for filtered VCFs (population structure, ps)
process qc_ps_post {

  input:
  tuple file(ps_vcf), file(ps_csi), file(gs_vcf), file(gs_csi)
  path windows_dir
  each window

  output:
  path "${ps_vcf.simpleName}_${window}_qc_ps_post.txt"

  script:
  """
  if [[ "${window}" == "scaffold"* ]];
  then
    bcftools stats --threads ${task.cpus} --regions-file ${windows_dir}/${window} ${ps_vcf} \
    > ${ps_vcf.simpleName}_${window}_qc_ps_post.txt
  else
    bcftools stats --threads ${task.cpus} --regions ${window} ${ps_vcf} \
    > ${ps_vcf.simpleName}_${window}_qc_ps_post.txt
  fi
  """
}

// Run QC for filtered VCFs (genome scan, gs)
process qc_gs_post {

  input:
  tuple file(ps_vcf), file(ps_csi), file(gs_vcf), file(gs_csi)
  path windows_dir
  each window

  output:
  path "${gs_vcf.simpleName}_${window}_qc_gs_post.txt"

  script:
  """
  if [[ "${window}" == "scaffold"* ]];
  then
    bcftools stats --threads ${task.cpus} --regions-file ${windows_dir}/${window} ${gs_vcf} \
    > ${gs_vcf.simpleName}_${window}_qc_gs_post.txt
  else
    bcftools stats --threads ${task.cpus} --regions ${window} ${gs_vcf} \
    > ${gs_vcf.simpleName}_${window}_qc_gs_post.txt
  fi
  """
}

// Combine output of QC processes
process combine_qc {

  publishDir "${params.publish_dir}/qc", saveAs: { filename -> "$filename" }, mode: 'copy'

  input:
  path qc_pre_files
  path qc_post_ps_files
  path qc_post_gs_files
  
  output:
  path "stats_pre.txt"
  path "stats_ps_post.txt"
  path "stats_gs_post.txt"

  script:
  """
  cat ${qc_pre_files} > stats_pre.txt
  cat ${qc_post_ps_files} > stats_ps_post.txt
  cat ${qc_post_gs_files} > stats_gs_post.txt
  """
}
