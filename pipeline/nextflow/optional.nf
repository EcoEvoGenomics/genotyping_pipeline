process deduplicate_reads {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory 4.GB
    time 1.h
    
    input:
    tuple val(sample), file(r1_reads), file(r2_reads), path(qcmetrics, stageAs: './qc-metrics/')

    output:
    tuple val(sample), file("${r1_reads.simpleName}.fastq.gz"), file("${r2_reads.simpleName}.fastq.gz"), path("qc-metrics/*")

    script:
    """
    seqkit rmdup -j ${task.cpus} --by-name -o ${r1_reads.simpleName}_deduplicated.fastq.gz ${r1_reads}
    seqkit rmdup -j ${task.cpus} --by-name -o ${r2_reads.simpleName}_deduplicated.fastq.gz ${r2_reads}
    mv ${r1_reads.simpleName}_deduplicated.fastq.gz ${r1_reads.simpleName}.fastq.gz
    mv ${r2_reads.simpleName}_deduplicated.fastq.gz ${r2_reads.simpleName}.fastq.gz
    seqkit stats -j ${task.cpus} -To qc-metrics/deduplication.tsv *.fastq.gz
    """
}

process downsample_reads {
    
    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory { 1.GB * Math.ceil(Math.max(r1_reads.size(), r2_reads.size()) / 1024 ** 3) }
    time 1.h

    input:
    tuple val(sample), file(r1_reads), file(r2_reads), path('qc-metrics/*', stageAs: './qc-metrics/')
    
    output:
    tuple val(sample), file("${r1_reads.simpleName}.fastq.gz"), file("${r2_reads.simpleName}.fastq.gz"), path("qc-metrics/*")

    script:
    """
    seqkit sample --two-pass -n ${params.read_target} -o ${r1_reads.simpleName}_downsampled.fastq.gz ${r1_reads}
    seqkit sample --two-pass -n ${params.read_target} -o ${r2_reads.simpleName}_downsampled.fastq.gz ${r2_reads}
    mv ${r1_reads.simpleName}_downsampled.fastq.gz ${r1_reads.simpleName}.fastq.gz
    mv ${r2_reads.simpleName}_downsampled.fastq.gz ${r2_reads.simpleName}.fastq.gz
    seqkit stats -j ${task.cpus} -To qc-metrics/downsampling.tsv *.fastq.gz
    """
}