#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Filter VCF
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default parameters
params.publish_dir = './output'
params.vcf_dir = "${params.publish_dir}/02-variants_unfiltered"

// Default filtering parameters
params.min_alleles=2
params.max_alleles=2
params.miss=0.8
params.q_site=30
params.q_site_gs=30
params.min_depth=5
params.max_depth=30
params.min_geno_depth=5
params.max_geno_depth=30
params.keep=""
params.stats_downsample_sites=10000

// Include duplicate processes
include { summarise_vcf; concatenate_all } from './call_variants.nf'

// Workflow
workflow{

  def filtered_chromosome_vcfs = Channel
  .fromPath("${params.vcf_dir}/**.vcf.gz")
  .map { vcf -> 
    def key = vcf.simpleName
    def index = vcf.toString().replace('vcf.gz', 'vcf.gz.csi')
    tuple(key, file(vcf), file(index))
  } \
  | filter_vcf

  def filtered_chromosome_vchks = filtered_chromosome_vcfs \
  | summarise_vcf

  concatenate_all(
        (filtered_chromosome_vcfs.flatten().collect()),
        (filtered_chromosome_vchks.collect()),
        "variants_${params.filtering_label}"
    )

}

// Filter a VCF
process filter_vcf {

  publishDir "${params.publish_dir}/chroms/${key}", saveAs: { filename -> "$filename" }, mode: 'copy'
  
  input:
  tuple val(key), path('input.vcf.gz'), path('input.vcf.gz.csi') 

  output:
  tuple \
  file("${key}_${params.filtering_label}.vcf.gz"), \
  file("${key}_${params.filtering_label}.vcf.gz.csi")

  script:
  """
  # FILTER WITH OR WITHOUT KEEP PARAMETER
  if [[ -f ${params.keep} ]]; then
    vcftools --gzvcf input.vcf.gz \
    --min-alleles ${params.min_alleles} \
    --max-alleles ${params.max_alleles} \
    --max-missing ${params.miss} \
    --minQ ${params.q_site} \
    --min-meanDP ${params.min_depth} \
    --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} \
    --maxDP ${params.max_geno_depth} \
    --keep ${params.keep} \
    --remove-filtered-all \
    --remove-indels \
    --recode-INFO-all \
    --recode \
    --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z \
    -o ${key}_${params.filtering_label}.vcf.gz
  else
    vcftools --gzvcf input.vcf.gz \
    --min-alleles ${params.min_alleles} \
    --max-alleles ${params.max_alleles} \
    --max-missing ${params.miss} \
    --minQ ${params.q_site} \
    --min-meanDP ${params.min_depth} \
    --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} \
    --maxDP ${params.max_geno_depth} \
    --remove-filtered-all \
    --remove-indels \
    --recode-INFO-all \
    --recode \
    --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z \
    -o ${key}_${params.filtering_label}.vcf.gz
  fi

  # INDEX FILTERED VCF
  bcftools index --threads ${task.cpus} ${key}_${params.filtering_label}.vcf.gz
  """
}
