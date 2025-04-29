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
        .multiMap { cols -> input_reads: [cols[0], cols[1], cols[2]] }
        .set { samples }

    def parsed_reads = parse_input_reads(samples.input_reads)

    if (params.deduplicate == 'yes') {
        parsed_reads = deduplicate_reads(parsed_reads)
    }

    if (params.downsample == 'yes') {
        parsed_reads = downsample_reads(parsed_reads)
    }

    def trimmed_reads = trim_reads(parsed_reads)
    def aligned_reads = align_reads(trimmed_reads, file(params.ref_genome), file(params.ref_index))

}

// Step 1 - Parse input reads for a sample
process parse_input_reads {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory 1.GB
    time 15.m

    input:
    tuple val(ID), path(R1), path(R2)

    output:
    tuple val(ID), path(R1), path(R2), path('qc-metrics/*')

    script:
    """
    mkdir qc-metrics
    seqkit stats -j ${task.cpus} -To qc-metrics/unprocessed.tsv *.fastq.gz
    printf '%s\\t%s\\n' \
        'ID' '${ID}' \
        'deduplicate' '${params.deduplicate}' \
        'downsample' '${params.downsample}' \
        > qc-metrics/settings.tsv
    """
}

// Step 2 - Read trim_reads
process trim_reads {

    publishDir "${params.publish_dir}/${ID}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    cpus 4
    memory 4.GB
    time 30.m

    input:
    tuple val(ID), file(R1), file(R2), path(qcmetrics, stageAs: './qc-metrics/')

    output:
    tuple val(ID), file("${ID}_R1.fastq.gz"), file("${ID}_R2.fastq.gz"), path("qc-metrics/*")

    script:
    """
    fastp \
    --in1 ${R1} \
    --in2 ${R2} \
    --out1 ${ID}_R1.fastq.gz \
    --out2 ${ID}_R2.fastq.gz \
    --report_title "${ID}" \
    --html qc-metrics/${ID}.html \
    --json qc-metrics/${ID}.json
    """
}

// Step 3 - Align to reference genome using GPU
process align_reads {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--nv"
    memory { 16.GB + 7.5.GB * Math.ceil(Math.max(r1_file.size(), r2_file.size()) / 1024 ** 3) }
    time 1.h
    
    label "gpu"

    input:
    tuple val(sample), path(r1_file), path(r2_file), path(qcmetrics, stageAs: './qc-metrics/')
    path('ref_genome.fa')
    path('ref_genome.fa.fai')

    output:
    tuple \
    val(sample), \
    path("${sample}.cram"), \
    path("${sample}.cram.crai"), \
    path('qc-metrics/*')

    script:
    """
    # Prepare QC directory
    mkdir qc-metrics

    # Align with Parabricks FQ2BAM
    pbrun fq2bam \
    --ref ref_genome.fa \
    --in-fq ${r1_file} ${r2_file} \
    --out-bam ${sample}.cram \
    --out-duplicate-metrics qc-metrics/dedup.txt \
    --out-qc-metrics-dir qc-metrics \
    --tmp .

    # Parse QC files
    rm qc-metrics/*.pdf
    rm qc-metrics/*.png
    for file in qc-metrics/*.txt; do
        mv \$file qc-metrics/${sample}.\$(basename \$file .txt)
    done
    """
}
}
