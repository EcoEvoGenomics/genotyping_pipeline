#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: MultiQC
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

workflow {

    def config = file("./pipeline/config/multiqc_config.yaml")
    def sparrows_logo = file("./pipeline/assets/sparrows.jpg")
    run_multiqc(params.results_dir, config, sparrows_logo)

}

process run_multiqc {

    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'
    
    input:
    path(results_dir)
    path(multiqc_config)
    path(sparrows_logo)

    output:
    file('genotyping_pipeline_multiqc_report.html')

    script:
    """
    mkdir multiqc_outputs
    multiqc \
    --config ${multiqc_config} \
    --force \
    ${results_dir}
    """
}
