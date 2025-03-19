#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Downsample BAMs
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default parameters
params.depth = 10
params.publish_dir = './output'

workflow{    
    bams_list = file(params.bams).readLines()
    Channel.fromList(bams_list)
        .set{bams}

    align_downsample(bams) | view
}

process align_downsample {

    publishDir "${params.publish_dir}/align_downsample", saveAs: { filename -> "$filename" }, mode: 'copy'
    
    input:
    path (bam)

    output:
    tuple stdout, path("${bam.simpleName}_ds.bam"), path("${bam.simpleName}_ds.bam.bai")

    script:
    """
    mean_coverage=\$(samtools depth ${bam} | awk '{sum += \$3} END {result = sum / NR; printf "%d\\n", result}' )
    fraction_sampled=\$(echo "scale=1; ${params.depth}/\${mean_coverage}" | bc)

    if (( \${mean_coverage} < ${params.depth} )); then
        echo "Depth for ${bam} is \${mean_coverage}, lower than ${params.depth}. Downsampled file will be identical to input."
    else
        echo "Downsampling ${bam} by \${fraction_sampled} ..."
    fi

    samtools view -bs \${fraction_sampled} ${bam} > ${bam.simpleName}_ds.cram
    samtools index ${bam.simpleName}_ds.cram
    """
}
