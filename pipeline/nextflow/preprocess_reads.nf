#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Preprocess reads
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

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

// Step 0 - Parse sample
process parse_sample {

    input:
    tuple val(sample), path(f_read), path(r_read)

    output:
    val(sample)
    path(f_read)
    path(r_read)
    path('qc-metrics/*')

    script:
    """
    mkdir qc-metrics
    cat "sample: ${sample}" > qc-metrics/${sample}.txt
    cat "deduplicate: ${params.deduplicate}" >> qc-metrics/${sample}.txt
    cat "downsample: ${params.downsample}" >> qc-metrics/${sample}.txt
    """
}

// Step 1 (Optional) - Deduplicate reads
process deduplicate_reads {
    
    input: 
    val(sample)
    file(r1_reads)
    file(r2_reads)
    file('qc-metrics/*')

    output:
    val(sample)
    file("${r1_reads.simpleName}_rmdup.fastq.gz")
    file("${r2_reads.simpleName}_rmdup.fastq.gz")
    file("qc-metrics/*")

    script:
    """
    seqkit stats -j ${task.cpus} -To qc-metrics/before_dedup.tsv *.fastq.gz
    seqkit rmdup -j ${task.cpus} --by-name -o ${r1_reads.simpleName}_rmdup.fastq.gz ${r1_reads}
    seqkit rmdup -j ${task.cpus} --by-name -o ${r2_reads.simpleName}_rmdup.fastq.gz ${r2_reads}
    seqkit stats -j ${task.cpus} -To qc-metrics/after_dedup.tsv *rmdup.fastq.gz
    """
}

// Step 2 (Optional) - Downsample reads
process downsample_reads {
    
    input: 
    val(sample)
    file(r1_reads)
    file(r2_reads)
    file('qc-metrics/*')
    
    output:
    val(sample)
    file("${r1_reads.simpleName}_sample.fastq.gz")
    file("${r2_reads.simpleName}_sample.fastq.gz")
    file("qc-metrics/*")

    script:
    """
    seqkit stats -j ${task.cpus} -To qc-metrics/before_downsample.tsv *.fastq.gz
    seqkit sample --two-pass -n ${params.read_target} -o ${r1_reads.simpleName}_sample.fastq.gz ${r1_reads}
    seqkit sample --two-pass -n ${params.read_target} -o ${r2_reads.simpleName}_sample.fastq.gz ${r2_reads}
    seqkit stats -j ${task.cpus} -To qc-metrics/after_downsample.tsv *sample.fastq.gz
    """
}

// Step 3 - Read trim_reads
process trim_reads {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input: 
    val(sample)
    file(r1_reads)
    file(r2_reads)
    file('qc-metrics/*')

    output:
    val(sample)
    file("${sample}_R1.fastq.gz")
    file("${sample}_R2.fastq.gz")
    file("qc-metrics/*")

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
