#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Downsample CRAMs
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default parameters
params.ref = file('/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa')
params.depth = 10
params.publish_dir = './output'

workflow{    
    crams_list = file(params.crams).readLines()
    Channel.fromList(crams_list)
        .set{crams}

    align_downsample(crams)
}

process align_downsample {

    publishDir "${params.publish_dir}/align_maxdepth", saveAs: { filename -> "$filename" }, mode: 'copy'
    
    input:
    path (cram)

    output:
    tuple path("${cram.simpleName}_maxdepth.cram"), path("${cram.simpleName}_maxdepth.cram.crai")

    script:
    """
    mean_coverage=\$(samtools depth ${cram} | awk '{sum += \$3} END {result = sum / NR; printf "%d\\n", result}' )
    fraction_sampled=\$(echo "scale=1; ${params.depth}/\${mean_coverage}" | bc)

    if (( \${mean_coverage} < ${params.depth} )); then
        echo "Depth for ${cram} is \${mean_coverage}, lower than ${params.depth}. Downsampled file will be identical to input."
    else
        echo "Downsampling ${cram} by \${fraction_sampled} ..."
    fi

    samtools view -@ ${task.cpus} -s \${fraction_sampled} -C -T ${params.ref} -o ${cram.simpleName}_maxdepth.cram ${cram}
    samtools index -@ ${task.cpus} ${cram.simpleName}_maxdepth.cram
    """
}
