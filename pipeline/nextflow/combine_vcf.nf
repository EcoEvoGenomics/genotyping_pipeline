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

// Workflow
workflow  {
   Channel.fromPath("${params.input_dir}/*${params.pop_structure_label}.vcf.gz").set{pop_structure}
   Channel.fromPath("${params.input_dir}/*${params.genome_scan_label}.vcf.gz").set{genome_scan}
   combine_vcf(pop_structure) | renormalise_vcf
   combine_vcf(genome_scan) | renormalise_vcf
}

process combine_vcf {

    publishDir "${params.publish_dir}/combined_vcf", saveAs: { filename -> "$filename" }

    input:
    path list_of_vcfs, stageAs: "staged/*"
    
    output:
    path "${list_of_vcfs.name()}.vcf.gz"

    script:
    """
    bcftools concat --threads 8 -n -O z -o ${list_of_vcfs.name()}.vcf.gz staged/*
    bcftools index ${list_of_vcfs.name()}.vcf.gz
    """
}

process renormalise_vcf {

    publishDir "${params.publish_dir}/combined_vcf", saveAs: { filename -> "$filename" }

    input:
    path vcf

    output:
    path "${vcf.baseName}_norm.vcf.gz"

    script:
    """
    bcftools norm -d none -O z -o ${vcf.baseName}_norm.vcf.gz ${vcf}
    bcftools index ${vcf.baseName}_norm.vcf.gz
    """
}
