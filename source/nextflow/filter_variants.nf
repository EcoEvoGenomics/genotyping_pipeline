include { summarise_vcf; concatenate_vchks; concatenate_vcfs } from './call_variants.nf'

workflow{

  def ref_index = file(params.ref_genome.toString() + ".fai")

  // Obtain chromosome-level unfiltered VCFs
  def chromosome_vcfs = Channel
  .fromPath("${params.vcf_dir}/**.vcf.gz")
  .map { vcf -> 
    def key = vcf.simpleName
    def index = vcf.toString().replace('vcf.gz', 'vcf.gz.csi')
    tuple(key, file(vcf), file(index))
  }

  // Filter, retain only individuals in keepfile
  def filtered_chromosome_vcfs = filter_vcf(chromosome_vcfs, file(params.filtering_flags))

  // Obtain summary stats chromosome-level VCF
  def filtered_chromosome_vchks = filtered_chromosome_vcfs \
  | summarise_vcf

  // Concatenate and output chromosome-level VCFs and VCHKs
  concatenate_vchks(filtered_chromosome_vchks.collect(), "${params.filtering_label}")
  concatenate_vcfs(filtered_chromosome_vcfs.flatten().collect(), ref_index, params.ref_scaffold_name, "${params.filtering_label}")

  // Separately:
  save_filters_to_file(file(params.filtering_flags))

}

// Filter a VCF
process filter_vcf {

  publishDir "${params.publish_dir}/chroms/${key}", saveAs: { filename -> "$filename" }, mode: 'copy'

  // Container: https://wave.seqera.io/view/builds/bd-39fc8ab24f49f2d6_1
  container "community.wave.seqera.io/library/bcftools_vcftools:39fc8ab24f49f2d6"
  cpus 2
  memory { 2.GB * task.attempt }
  time { 4.h * task.attempt }
  
  errorStrategy "retry"
  maxRetries 3
  
  input:
  tuple val(key), path('input.vcf.gz'), path('input.vcf.gz.csi')
  path(filterfile)

  output:
  tuple \
  path("${key}_${params.filtering_label}.vcf.gz"), \
  path("${key}_${params.filtering_label}.vcf.gz.csi")

  script:
  """
  echo '--gzvcf input.vcf.gz' >> vcftools.args
  echo '--remove-indels' >> vcftools.args
  echo '--recode-INFO-all' >> vcftools.args
  echo '--recode' >> vcftools.args
  echo '--stdout' >> vcftools.args
  cat ${filterfile} >> vcftools.args

  cat vcftools.args \
  | xargs vcftools \
  | bcftools view --threads ${task.cpus} \
    -e 'ALT="*" || TYPE!="snp"' \
    -O z -o ${key}_${params.filtering_label}.vcf.gz

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

  input:
  path(filters)

  output:
  path("vcftools_${params.filtering_label}.txt")

  script:
  """
  mv ${filters} vcftools_${params.filtering_label}.txt
  """
}
