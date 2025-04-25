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
    seqkit stats -j ${task.cpus} -To qc-metrics/unprocessed.tsv *.fastq.gz
    printf '%s\\t%s\\n' \
        'sample' '${sample}' \
        'deduplicate' '${params.deduplicate}' \
        'downsample' '${params.downsample}' \
        > qc-metrics/settings.tsv
    """
}

// Step 1 (Optional) - Deduplicate reads
process deduplicate_reads {
    
    input: 
    val(sample)
    file(r1_reads)
    file(r2_reads)
    path('qc-metrics/*'), stageAs: './qc-metrics/'

    output:
    val(sample)
    file("${r1_reads.simpleName}.fastq.gz")
    file("${r2_reads.simpleName}.fastq.gz")
    path("qc-metrics/*")

    script:
    """
    seqkit rmdup -j ${task.cpus} --by-name -o ${r1_reads.simpleName}_deduplicated.fastq.gz ${r1_reads}
    seqkit rmdup -j ${task.cpus} --by-name -o ${r2_reads.simpleName}_deduplicated.fastq.gz ${r2_reads}
    mv ${r1_reads.simpleName}_deduplicated.fastq.gz ${r1_reads.simpleName}.fastq.gz
    mv ${r2_reads.simpleName}_deduplicated.fastq.gz ${r2_reads.simpleName}.fastq.gz
    seqkit stats -j ${task.cpus} -To qc-metrics/deduplication.tsv *.fastq.gz
    """
}

// Step 2 (Optional) - Downsample reads
process downsample_reads {
    
    memory { 1.GB * Math.ceil(Math.max(r1_reads.size(), r2_reads.size()) / 1024 ** 3) }

    input: 
    val(sample)
    file(r1_reads)
    file(r2_reads)
    path('qc-metrics/*'), stageAs: './qc-metrics/'
    
    output:
    val(sample)
    file("${r1_reads.simpleName}.fastq.gz")
    file("${r2_reads.simpleName}.fastq.gz")
    path("qc-metrics/*")

    script:
    """
    seqkit sample --two-pass -n ${params.read_target} -o ${r1_reads.simpleName}_downsampled.fastq.gz ${r1_reads}
    seqkit sample --two-pass -n ${params.read_target} -o ${r2_reads.simpleName}_downsampled.fastq.gz ${r2_reads}
    mv ${r1_reads.simpleName}_downsampled.fastq.gz ${r1_reads.simpleName}.fastq.gz
    mv ${r2_reads.simpleName}_downsampled.fastq.gz ${r2_reads.simpleName}.fastq.gz
    seqkit stats -j ${task.cpus} -To qc-metrics/downsampling.tsv *.fastq.gz
    """
}

// Step 3 - Read trim_reads
process trim_reads {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input: 
    val(sample)
    file("${sample}_R1.fastq.gz")
    file("${sample}_R2.fastq.gz")
    path('qc-metrics/*'), stageAs: './qc-metrics/'

    output:
    val(sample)
    file("${sample}_R1_TRIM.fastq.gz")
    file("${sample}_R2_TRIM.fastq.gz")
    path("qc-metrics/*")

    script:
    """
    fastp \
    --in1 ${sample}_R1.fastq.gz \
    --in2 ${sample}_R2.fastq.gz \
    --out1 ${sample}_R1_TRIM.fastq.gz \
    --out2 ${sample}_R2_TRIM.fastq.gz \
    --report_title "${sample}" \
    --html qc-metrics/${sample}.html \
    --json qc-metrics/${sample}.json
    """
}
