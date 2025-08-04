#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Filter VCF
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Include duplicate processes
include { summarise_vcf; concatenate_vchks; concatenate_vcfs } from './call_variants.nf'

// Workflow
workflow{

  // Obtain chromosome-level unfiltered VCFs and filter
  def filtered_chromosome_vcfs = Channel
  .fromPath("${params.vcf_dir}/**.vcf.gz")
  .map { vcf -> 
    def key = vcf.simpleName
    def index = vcf.toString().replace('vcf.gz', 'vcf.gz.csi')
    tuple(key, file(vcf), file(index))
  } \
  | filter_vcf

  // Obtain summary stats chromosome-level VCF
  def filtered_chromosome_vchks = filtered_chromosome_vcfs \
  | summarise_vcf

  // Concatenate and output chromosome-level VCFs and VCHKs
  concatenate_vchks(filtered_chromosome_vchks.collect(), "variants_${params.filtering_label}")
  concatenate_vcfs(filtered_chromosome_vcfs.flatten().collect(), params.ref_index, params.ref_scaffold_name, "variants_${params.filtering_label}")

  // Separately:
  save_filters_to_file()

}

// Filter a VCF
process filter_vcf {

  publishDir "${params.publish_dir}/chroms/${key}", saveAs: { filename -> "$filename" }, mode: 'copy'

  // Container: https://wave.seqera.io/view/builds/bd-39fc8ab24f49f2d6_1
  container "community.wave.seqera.io/library/bcftools_vcftools:39fc8ab24f49f2d6"
  cpus 2
  memory 1.GB
  time 2.h
  
  input:
  tuple val(key), path('input.vcf.gz'), path('input.vcf.gz.csi') 

  output:
  tuple \
  path("${key}_${params.filtering_label}.vcf.gz"), \
  path("${key}_${params.filtering_label}.vcf.gz.csi")

  script:
  """
  # FILTER WITH OR WITHOUT KEEP PARAMETER
  if [[ -f ${params.keep} ]]; then
    vcftools --gzvcf input.vcf.gz \
    --min-alleles ${params.min_alleles} \
    --max-alleles ${params.max_alleles} \
    --max-missing ${params.max_missing} \
    --min-meanDP ${params.min_meanDP} \
    --max-meanDP ${params.max_meanDP} \
    --minDP ${params.minDP} \
    --maxDP ${params.maxDP} \
    --minQ ${params.minQ} \
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
    --max-missing ${params.max_missing} \
    --min-meanDP ${params.min_meanDP} \
    --max-meanDP ${params.max_meanDP} \
    --minDP ${params.minDP} \
    --maxDP ${params.maxDP} \
    --minQ ${params.minQ} \
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

// Save the filters applied to a file
process save_filters_to_file {

  publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'

  cpus 1
  memory 256.MB
  time 5.m

  output:
  path("vcftools_${params.filtering_label}.tsv")

  script:
  """
  printf '%s\\t%s\\n' \
    'min-alleles' '${params.min_alleles}' \
    'max-alleles' '${params.max_alleles}' \
    'max-missing' '${params.max_missing}' \
    'min-meanDP' '${params.min_meanDP}' \
    'max-meanDP' '${params.max_meanDP}' \
    'minDP' '${params.minDP}' \
    'maxDP' '${params.maxDP}' \
    'minQ' '${params.minQ}' \
    'keep' '${params.keep}' \
    > vcftools_${params.filtering_label}.tsv
  """
}
