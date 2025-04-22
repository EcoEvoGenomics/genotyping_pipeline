#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: MultiQC
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

workflow {

    def config = file(params.multiqc_config)
    run_multiqc(params.results_dir, config)

}

process run_multiqc {

    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'
    
    input:
    path(results_dir)
    path(multiqc_config)

    output:
    file('multiqc_report.html')

    script:
    """
    mkdir multiqc_outputs
    multiqc \
    --config ${multiqc_config} \
    --force \
    ${results_dir}
    """
}
