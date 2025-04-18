#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Trim reads
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

workflow {
    
    // Input 'samples' is .CSV with columns for ID, F_READ_PATH, R_READ_PATH
    Channel.fromPath(params.samples)
        .splitCsv()
        .multiMap { cols -> sample: [cols[0], cols[1], cols[2]] }
        .set { samples }

    parse_sample(samples.sample) | trim_reads

}

// Step 0 - Parse sample
process parse_sample {

    input:
    tuple val(sample), path(f_read), path(r_read)

    output:
    val(sample)
    path(f_read)
    path(r_read)

    script:
    """
    echo "Initiating pipeline for ${sample} ..."
    """
}

// Step 1 - Read trim_reads
process trim_reads {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input: 
    val(sample)
    file("${sample}_R1.fastq.gz")
    file("${sample}_R2.fastq.gz")

    output:
    val(sample)
    file("${sample}_R1_TRIM.fastq.gz")
    file("${sample}_R2_TRIM.fastq.gz")
    file("${sample}.html")
    file("${sample}.json")

    script:
    """
    fastp \
    --in1 ${sample}_R1.fastq.gz \
    --in2 ${sample}_R2.fastq.gz \
    --out1 ${sample}_R1_TRIM.fastq.gz \
    --out2 ${sample}_R2_TRIM.fastq.gz \
    --report_title "${sample}" \
    --html ${sample}.html \
    --json ${sample}.json
    """
}
