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
params.vcf_dir = "${params.publish_dir}/called_variants"

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

  def vcf_and_index = Channel
  .fromPath("${params.vcf_dir}/vcf/*.vcf.gz")
  .map { vcf -> 
    def index = vcf.toString().replace('vcf.gz', 'vcf.gz.csi')
    tuple(file(vcf), file(index))
  }
 
  nonfilt(vcf_and_index)
  popfilt(vcf_and_index)
  genfilt(vcf_and_index)

}

workflow nonfilt {
  take: vcf_and_index
  main: vcf_and_index | downsample_vcf | get_vcf_stats | plot_vcf_stats
}

workflow popfilt {
  take: vcf_and_index
  main: vcf_and_index | vcf_filter_pop_structure | downsample_vcf | get_vcf_stats | plot_vcf_stats
}

workflow genfilt {
  take: vcf_and_index
  main: vcf_and_index | vcf_filter_genome_scan | downsample_vcf | get_vcf_stats | plot_vcf_stats
}

// Filtering for population structure analyses ...
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

// Randomly subsample a VCF for less intensive stats calculations
process downsample_vcf {

  input:
  tuple file(vcf), file(csi)

  output:
  tuple file("${vcf.simpleName}_downsampled.vcf.gz"), file("${vcf.simpleName}_downsampled.vcf.gz.csi")

  script:
  """
  # Ensure number of sampled sites does not exceed sites in VCF (currently redundant)
  sampled_sites=${params.stats_downsample_sites}
  vcf_num_sites=\$(bcftools view -H ${vcf} | wc -l)
  if (( \$sampled_sites > \$vcf_num_sites )); then
    sampled_sites=\$vcf_num_sites
  fi
  
  # Then, downsample ...
  bcftools view -Ov ${vcf} | vcfrandomsample -r 0.95 > ${vcf.simpleName}_downsampled.vcf
  bgzip ${vcf.simpleName}_downsampled.vcf
  bcftools index ${vcf.simpleName}_downsampled.vcf.gz
  """
}

// Calculate VCF stats
process get_vcf_stats {

  publishDir "${params.publish_dir}/stats", saveAs: { filename -> "$filename" }, mode: 'copy'

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
process plot_vcf_stats {
  
  publishDir "${params.publish_dir}/plots", saveAs: { filename -> "$filename" }, mode: 'copy'

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
  file "${vcf_name}/*"

  script:
  """
  mkdir ${vcf_name}

  Rscript -e \"
  library(ggplot2)

  var_qual <- read.table('${vcf_stats_lqual}', sep = '\\t', col.names = c('chr', 'pos', 'qual'), header = TRUE)
  a <- ggplot(var_qual, aes(qual)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/variant_quality.png', plot = a + theme_light())

  var_depth <- read.table('${vcf_stats_ldepth_mean}', sep = '\\t', col.names = c('chr', 'pos', 'mean_depth', 'var_depth'), header = TRUE)
  a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/variant_mean_depth.png', plot = a + theme_light())

  var_miss <- read.table('${vcf_stats_lmiss}', sep = '\\t', col.names = c('chr', 'pos', 'nchr', 'nfiltered', 'nmiss', 'fmiss'), header = TRUE)
  a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/variant_missingness.png', plot = a + theme_light())

  var_freq <- read.table('${vcf_stats_frq}', sep = '\\t', col.names = c('chr', 'pos', 'nalleles', 'nchr', 'a1', 'a2'), fill = TRUE, skip = 1)
  var_freq[["maf"]] <- apply(var_freq[c('a1', 'a2')], 1, function(z) min(z))
  a <- ggplot(var_freq, aes(maf)) + geom_density(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/minor_allele_frequency.png', plot = a + theme_light())
  
  ind_depth <- read.table('${vcf_stats_idepth}', sep = '\\t', col.names = c('ind', 'nsites', 'depth'), header = TRUE)
  a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/individual_depth.png', plot = a + theme_light())

  ind_miss  <- read.table('${vcf_stats_imiss}', sep = '\\t', col.names = c('ind', 'ndata', 'nfiltered', 'nmiss', 'fmiss'), header = TRUE)
  a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/individual_missingness.png', plot = a + theme_light())

  ind_het <- read.table('${vcf_stats_het}', sep = '\\t', col.names = c('ind','ho', 'he', 'nsites', 'f'), header = TRUE)
  a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = 'dodgerblue1', colour = 'black', alpha = 0.3)
  ggsave('${vcf_name}/individual_het.png', plot = a + theme_light())
  \"
  """
}
