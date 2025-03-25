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
params.vcf_dir = "${params.publish_dir}/raw_vcf"

// Default filtering parameters
params.miss=0.8
params.q_site_ps=30
params.q_site_gs=30
params.min_depth=5
params.max_depth=30
params.min_geno_depth=5
params.max_geno_depth=30
params.keep=""
params.stats_downsample_sites=10000

// Workflow
workflow{

  def raw_vcf_and_index = Channel
  .fromPath("${params.vcf_dir}/vcf/*.vcf.gz")
  .map { vcf -> 
    def index = vcf.toString().replace('vcf.gz', 'vcf.gz.csi')
    tuple(file(vcf), file(index))
  }

  raw_vcf_and_index | downsample_vcf | get_vcf_stats | analyse_vcf_stats

  def indels_removed = raw_vcf_and_index | rm_indels
  filter_pop_structure(indels_removed)
  filter_genome_structure(indels_removed)

}

workflow filter_pop_structure {
  
  take:
  indels_removed

  main:
  indels_removed | vcf_filter_pop_structure | downsample_vcf | get_vcf_stats | analyse_vcf_stats

}

workflow filter_genome_structure {
  
  take:
  indels_removed

  main:
  indels_removed | vcf_filter_genome_scan | downsample_vcf | get_vcf_stats | analyse_vcf_stats

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
  tuple file(vcf), file(vcf_index) 

  output:
  tuple file("${vcf.simpleName}_filt_ps.vcf.gz"), file("${vcf.simpleName}_filt_ps.vcf.gz.csi")

  script:
  """
  if [[ -f ${params.keep} ]]; then
    vcftools --gzvcf $vcf --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing ${params.miss} \
    --minQ ${params.q_site_ps} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${vcf.simpleName}_filt_ps.vcf.gz
  else
    vcftools --gzvcf $vcf --remove-indels --remove-filtered-all \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing ${params.miss} \
    --minQ ${params.q_site_ps} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${vcf.simpleName}_filt_ps.vcf.gz
  fi
  bcftools index --threads ${task.cpus} ${vcf.simpleName}_filt_ps.vcf.gz
  """
}

// ... and for genome scans.
process vcf_filter_genome_scan {

  publishDir "${params.publish_dir}/vcf_filtered", saveAs: { filename -> "$filename" }, mode: 'copy'

  input:
  tuple file(vcf), file(vcf_index) 

  output:
  tuple file("${vcf.simpleName}_filt_gs.vcf.gz"), file("${vcf.simpleName}_filt_gs.vcf.gz.csi")

  script:
  """
  if [[ -f ${params.keep} ]]; then
    vcftools --gzvcf $vcf --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --max-alleles 2 \
    --minQ ${params.q_site_gs} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${vcf.simpleName}_filt_gs.vcf.gz
  else
    vcftools --gzvcf $vcf --remove-indels --remove-filtered-all \
    --max-alleles 2 \
    --minQ ${params.q_site_gs} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout \
  | bcftools view --threads ${task.cpus} -e 'N_ALT>1' -O z -o ${vcf.simpleName}_filt_gs.vcf.gz
  fi
  bcftools index --threads ${task.cpus} ${vcf.simpleName}_filt_gs.vcf.gz
  """
}

// Randomly subsample a VCF
process downsample_vcf {

  input:
  tuple file(vcf), file(csi)

  output:
  tuple file("${vcf.simpleName}_downsampled.vcf.gz"), file("${vcf.simpleName}_downsampled.vcf.gz.csi")

  script:
  """
  # Ensure number of sampled sites does not exceed sites in VCF
  sampled_sites=${params.stats_downsample_sites}
  vcf_num_sites=\$(bcftools view -H ${vcf} | wc -l)
  if (( \$sampled_sites > \$vcf_num_sites )); then
    sampled_sites=\$vcf_num_sites
  fi
  
  # Then, downsample ...
  bcftools view --header-only ${vcf} > ${vcf.simpleName}_downsampled.vcf
  bcftools view --no-header ${vcf} \
    | awk '{printf("%f\\t%s\\n",rand(),\$0);}' \
    | sort -t \$'\\t'  -T . -k1, \
    | head -n \$sampled_sites \
    | cut -f 2- \
    >> ${vcf.simpleName}_downsampled.vcf
  bgzip ${vcf.simpleName}_downsampled.vcf
  bcftools index ${vcf.simpleName}_downsampled.vcf.gz
  """
}

// Calculate VCF stats
process get_vcf_stats {

  input:
  tuple path(vcf), path(csi)

  output:
  file "${vcf.simpleName}_stats.frq"
  file "${vcf.simpleName}_stats.idepth"
  file "${vcf.simpleName}_stats.ldepth.mean"
  file "${vcf.simpleName}_stats.lqual"
  file "${vcf.simpleName}_stats.imiss"
  file "${vcf.simpleName}_stats.lmiss"
  file "${vcf.simpleName}_stats.het"
  val vcf.simpleName

  script:
  """
  vcftools --gzvcf ${vcf} --freq2 --out ${vcf.simpleName}_stats --max-alleles 2
  vcftools --gzvcf ${vcf} --depth --out ${vcf.simpleName}_stats
  vcftools --gzvcf ${vcf} --site-mean-depth --out ${vcf.simpleName}_stats
  vcftools --gzvcf ${vcf} --site-quality --out ${vcf.simpleName}_stats
  vcftools --gzvcf ${vcf} --missing-indv --out ${vcf.simpleName}_stats
  vcftools --gzvcf ${vcf} --missing-site --out ${vcf.simpleName}_stats
  vcftools --gzvcf ${vcf} --het --out ${vcf.simpleName}_stats
  """
}

// Collate different VCF stats, plot, and output
process analyse_vcf_stats {
  
  publishDir "${params.publish_dir}/stats", saveAs: { filename -> "$filename" }, mode: 'copy'

  input:
  file vcf_stats_frq
  file vcf_stats_idepth
  file vcf_stats_ldepth_mean
  file vcf_stats_lqual
  file vcf_stats_imiss
  file vcf_stats_lmiss
  file vcf_stats_het
  val vcf_name

  output:
  file "analysis_plots/*"
  file vcf_stats_frq
  file vcf_stats_idepth
  file vcf_stats_ldepth_mean
  file vcf_stats_lqual
  file vcf_stats_imiss
  file vcf_stats_lmiss
  file vcf_stats_het

  script:
  """
  Rscript -e \"
  library(tidyverse)

  # Rename files to simplify processing in R
  file.copy('${vcf_stats_frq}', './vcf.frq')
  file.copy('${vcf_stats_idepth}', './vcf.idepth')
  file.copy('${vcf_stats_ldepth_mean}', './vcf.ldepth.mean')
  file.copy('${vcf_stats_lqual}', './vcf.lqual')
  file.copy('${vcf_stats_imiss}', './vcf.imiss')
  file.copy('${vcf_stats_lmiss}', './vcf.lmiss')
  file.copy('${vcf_stats_het}', './vcf.het')

  var_qual <- read_delim('./vcf.lqual', delim = '\\t', col_names = c('chr', 'pos', 'qual'), skip = 1)
  a <- ggplot(var_qual, aes(qual)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_variant_quality.png', plot = a + theme_light())

  var_depth <- read_delim('./vcf.ldepth.mean', delim = '\\t', col_names = c('chr', 'pos', 'mean_depth', 'var_depth'), skip = 1)
  a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_variant_mean_depth.png', plot = a + theme_light())

  var_miss <- read_delim('./vcf.lmiss', delim = '\\t', col_names = c('chr', 'pos', 'nchr', 'nfiltered', 'nmiss', 'fmiss'), skip = 1)
  a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_variant_missingness.png', plot = a + theme_light())

  var_freq <- read_delim('./vcf.frq', delim = '\\t', col_names = c('chr', 'pos', 'nalleles', 'nchr', 'a1', 'a2'), skip = 1)
  var_freq\$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
  a <- ggplot(var_freq,(maf)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_minor_allele_frequency.png', plot = a + theme_light())
  
  ind_depth <- read_delim('./vcf.idepth', delim = '\\t', col_names = c('ind', 'nsites', 'depth'), skip = 1)
  a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_individual_depth.png', plot = a + theme_light())

  ind_miss  <- read_delim('./vcf.imiss', delim = '\\t', col_names = c('ind', 'ndata', 'nfiltered', 'nmiss', 'fmiss'), skip = 1)
  a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_individual_missingness.png', plot = a + theme_light())

  ind_het <- read_delim('./vcf.het', delim = '\\t', col_names = c('ind','ho', 'he', 'nsites', 'f'), skip = 1)
  a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('analysis_plots/${vcf_name}_individual_het.png', plot = a + theme_light())
  \"
  """
}
