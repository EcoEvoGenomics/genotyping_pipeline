#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Downsample BAMs
//
// Developed by Mark Ravinet and Erik Sandertun Røed
// Maintained by Erik Sandertun Røed

// Default parameters
params.input_dir = '.output/filt_vcf/vcf_filtered'
params.pop_structure_label = '_ps'
params.genome_scan_label = '_gs'
params.publish_dir = './output/'

workflow  {
    population_structure()
    genome_scan()
}

workflow population_structure {
    Channel.fromPath("${params.input_dir}/*${params.pop_structure_label}.vcf.gz").collect().set{pop_structure}
    combine_vcf('combined_ps', pop_structure) | renormalise_vcf
}

workflow genome_scan {
    Channel.fromPath("${params.input_dir}/*${params.genome_scan_label}.vcf.gz").collect().set{genome_scan}
    combine_vcf('combined_gs', genome_scan) | renormalise_vcf
}

process combine_vcf {

    publishDir "${params.publish_dir}/combined_vcf", saveAs: { filename -> "$filename" }, mode: 'copy'
    input:
    val name
    path list_of_vcfs, stageAs: "staged/*"
    
    output:
    path "${name}.vcf.gz"

    script:
    """
    bcftools concat --threads 8 -n -O z -o ${name}.vcf.gz staged/*
    bcftools index ${name}.vcf.gz
    """
}

process renormalise_vcf {

    publishDir "${params.publish_dir}/combined_vcf", saveAs: { filename -> "$filename" }, mode: 'copy'
    input:
    path vcf

    output:
    path "${vcf.getName().replaceAll(/\.vcf\.gz$/, '')}_norm.vcf.gz"

    script:
    """
    vcfName="${vcf.getName().replaceAll(/\.vcf\.gz$/, '')}"
    bcftools norm -d none -O z -o \${vcfName}_norm.vcf.gz ${vcf}
    bcftools index \${vcfName}_norm.vcf.gz
    """
}
