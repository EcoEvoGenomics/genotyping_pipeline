#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: MultiQC
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

workflow {

    run_multiqc(params.results_dir, params.multiqc_config)

}

process run_multiqc {

    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'
    
    input:
    path(results_dir)
    file(multiqc_config)

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
