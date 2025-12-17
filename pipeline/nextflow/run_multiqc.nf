#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: MultiQC
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

workflow {

    def config = file("./pipeline/assets/multiqc_config.yaml")
    def sparrows_logo = file("./pipeline/assets/sparrows.jpg")
    run_multiqc(params.results_dir, config, sparrows_logo)

}

process run_multiqc {

    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'
    
    container "quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0"
    cpus 1
    memory 3.GB
    time 30.m

    input:
    path(results_dir)
    path(multiqc_config)
    path(sparrows_logo)

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc \
    --config ${multiqc_config} \
    --force \
    ${results_dir}
    mv *multiqc_report.html multiqc_report.html
    """
}
