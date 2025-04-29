#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Preprocess reads
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

include { downsample_reads; deduplicate_reads } from "./optional.nf"

workflow {
    
    // Input 'samples' is .CSV with columns for ID, F_READ_PATH, R_READ_PATH
    Channel.fromPath(params.samples)
        .splitCsv()
        .multiMap { cols -> sample: [cols[0], cols[1], cols[2]] }
        .set { samples }

    def parsed_samples = parse_sample(samples.sample)

    if (params.deduplicate == 'yes') {
        parsed_samples = parsed_samples | deduplicate_reads
    }

    if (params.downsample == 'yes') {
        parsed_samples = parsed_samples | downsample_reads
    }

    parsed_samples | trim_reads

}

// Step 1 - Parse sample
process parse_sample {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory 1.GB
    time 15.m

    input:
    tuple val(sample), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample), path(r1_reads), path(r2_reads), path('qc-metrics/*')

    script:
    """
    mkdir qc-metrics
    seqkit stats -j ${task.cpus} -To qc-metrics/unprocessed.tsv *.fastq.gz
    printf '%s\\t%s\\n' \
        'sample' '${sample}' \
        'deduplicate' '${params.deduplicate}' \
        'downsample' '${params.downsample}' \
        > qc-metrics/settings.tsv
    """
}

// Step 2 - Read trim_reads
process trim_reads {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    cpus 4
    memory 4.GB
    time 30.m

    input:
    tuple val(sample), file(r1_reads), file(r2_reads), path(qcmetrics, stageAs: './qc-metrics/')

    output:
    tuple val(sample), file("${sample}_R1.fastq.gz"), file("${sample}_R2.fastq.gz"), path("qc-metrics/*")

    script:
    """
    fastp \
    --in1 ${r1_reads} \
    --in2 ${r2_reads} \
    --out1 ${sample}_R1.fastq.gz \
    --out2 ${sample}_R2.fastq.gz \
    --report_title "${sample}" \
    --html qc-metrics/${sample}.html \
    --json qc-metrics/${sample}.json
    """
}
