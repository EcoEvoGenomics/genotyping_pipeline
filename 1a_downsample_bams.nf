#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Downsample BAMs
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// set depth cut off
params.depth = 10

// create bams channel
bams_list = file(params.bams)
    .readLines()
// create windows channel
Channel
    .fromList( bams_list )
    .set{bams}

//bams.view()

// Step 1 - calculate  statistics
process align_downsample {

    publishDir 'align_downsample', saveAs: { filename -> "$filename" }, mode: 'copy'
    
    input:
    path (bam)

    output:
    // stdout
    tuple stdout, path("${bam.simpleName}_ds.bam"), path("${bam.simpleName}_ds.bam.bai"), optional: true

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

// workflow starts here!

workflow{    
    align_downsample(bams) | view
}
