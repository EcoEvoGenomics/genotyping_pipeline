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
    tuple stdout, path("${bam.simpleName}_ds.bam"), path("${bam.simpleName}_ds.bam.bai"), optional: true

    script:
    """
    # work out mean coverage
    ##COV=\$(samtools depth ${bam} | awk '{sum += \$3} END {result = sum / NR; printf "%.2f\\n", result}' )
    COV=\$(samtools depth ${bam} | awk '{sum += \$3} END {result = sum / NR; printf "%d\\n", result}' )

    if (( \${COV} > ${params.depth} )); then

        echo "Depth for ${bam} is \${COV} - i.e. greater than ${params.depth} - downsampling."
        #FRAC=\$(expr ${params.depth} \* \${COV})
        FRAC=\$(echo "scale=1; ${params.depth}/\${COV}" | bc)
        echo "Will downsampled by \${FRAC}"
        samtools view -bs \${FRAC} ${bam} > ${bam.simpleName}_ds.bam
        samtools index ${bam.simpleName}_ds.bam

    else
        echo "Depth is \${COV} - i.e. less than ${params.depth} - no need to downsample."
    fi
    """
}
