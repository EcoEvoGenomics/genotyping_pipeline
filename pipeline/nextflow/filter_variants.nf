#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Filter VCF
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default parameters
params.vcf_dir = "${params.publish_dir}/raw_vcf"
params.publish_dir = './output'

// Default filtering parameters
params.miss=0.8
params.q_site_filt_ps=30
params.q_site_filt_gs=30
params.min_depth=5
params.max_depth=30
params.min_geno_depth=5
params.max_geno_depth=30
params.keep=""

// Workflow
workflow{

  def raw_vcf_and_index = Channel
  .fromPath("${params.vcf_dir}/vcf/*.vcf.gz")
  .map { vcf -> 
    def index = vcf.toString().replace('vcf.gz', 'vcf.gz.csi')
    tuple(file(vcf), file(index))
  }

  def indels_removed = raw_vcf_and_index | rm_indels
  def stats_post_filtering_pop_structure = indels_removed | vcf_filter_pop_structure | vcf_stats
  def stats_post_filtering_genome_scan = indels_removed | vcf_filter_genome_scan | vcf_stats
  def stats_pre_filtering = raw_vcf_and_index | vcf_stats

  collate_filtering_report(
    stats_pre_filtering,
    stats_post_filtering_pop_structure,
    stats_post_filtering_genome_scan
  )

}

// Filtering, Step 1 - Remove spanning indels and re-normalise
process rm_indels {

  input:
  tuple file(raw_vcf), file(raw_vcf_index)

  output:
  tuple \
    file ("${raw_vcf.simpleName}.vcf.gz"), \
    file ("${raw_vcf.simpleName}.vcf.gz.csi") 

  script:
  """
  bcftools view --threads ${task.cpus} -V indels -e 'ALT="*" | N_ALT>1' $raw_vcf \
    | bcftools norm --threads ${task.cpus} -D -O z -o ${raw_vcf.simpleName}.vcf.gz
  bcftools index --threads ${task.cpus} ${raw_vcf.simpleName}.vcf.gz
  """
}

// Filtering, Step 2 - Filtering for population structure analyses ...
process vcf_filter_pop_structure {
  publishDir "${params.publish_dir}/vcf_filtered", saveAs: { filename -> "$filename" }, mode: 'copy'
  
  input:
  tuple file(rm_indel_vcf), file(rm_indel_vcf_index) 

  output:
  tuple file("${rm_indel_vcf.simpleName}_filt_ps.vcf.gz"), file("${rm_indel_vcf.simpleName}_filt_ps.vcf.gz.csi")

  script:
  """
  if [[ -f ${params.keep} ]]; then
    vcftools --gzvcf $rm_indel_vcf --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing ${params.miss} \
    --minQ ${params.q_site_filt_ps} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filt_ps.vcf.gz
  else
    vcftools --gzvcf $rm_indel_vcf --remove-indels --remove-filtered-all \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing ${params.miss} \
    --minQ ${params.q_site_filt_ps} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filt_ps.vcf.gz
  fi
  bcftools index --threads ${task.cpus} ${rm_indel_vcf.simpleName}_filt_ps.vcf.gz
  """
}

// ... and for genome scans.
process vcf_filter_genome_scan {
  publishDir "${params.publish_dir}/vcf_filtered", saveAs: { filename -> "$filename" }, mode: 'copy'

  input:
  tuple file(rm_indel_vcf), file(rm_indel_vcf_index) 

  output:
  tuple file("${rm_indel_vcf.simpleName}_filt_gs.vcf.gz"), file("${rm_indel_vcf.simpleName}_filt_gs.vcf.gz.csi")

  script:
  """
  if [[ -f ${params.keep} ]]; then
    vcftools --gzvcf $rm_indel_vcf --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --max-alleles 2 \
    --minQ ${params.q_site_filt_gs} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filt_gs.vcf.gz
  else
    vcftools --gzvcf $rm_indel_vcf --remove-indels --remove-filtered-all \
    --max-alleles 2 \
    --minQ ${params.q_site_filt_gs} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${rm_indel_vcf.simpleName}_filt_gs.vcf.gz
  fi
  bcftools index --threads ${task.cpus} ${rm_indel_vcf.simpleName}_filt_gs.vcf.gz
  """
}

// Randomly subsample a VCF and calculate stats
process vcf_stats {

  input:
  tuple file(vcf), file(csi)

  output:
  path "${vcf.simpleName}_stats.txt"

  script:
  """
  
  """
}

process collate_filtering_report {

  input:
  file stats_pre_filtering
  file stats_post_filtering_pop_structure
  file stats_post_filtering_genome_scan
  
  output:
  file ('report.txt')

  script:
  """
  """
}
