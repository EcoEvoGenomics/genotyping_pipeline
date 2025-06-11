process deduplicate_reads {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory 4.GB
    time 1.h
    
    input:
    tuple val(ID), val(LANE), file(R1), file(R2), path(qcmetrics, stageAs: './qc-metrics/')

    output:
    tuple val(ID), val(LANE), file("${R1.simpleName}.fastq.gz"), file("${R2.simpleName}.fastq.gz"), path("qc-metrics/*")

    script:
    """
    seqkit rmdup -j ${task.cpus} --by-name -o ${R1.simpleName}_deduplicated.fastq.gz ${R1}
    seqkit rmdup -j ${task.cpus} --by-name -o ${R2.simpleName}_deduplicated.fastq.gz ${R2}
    mv ${R1.simpleName}_deduplicated.fastq.gz ${R1.simpleName}.fastq.gz
    mv ${R2.simpleName}_deduplicated.fastq.gz ${R2.simpleName}.fastq.gz
    seqkit stats -j ${task.cpus} -To qc-metrics/deduplication.tsv *.fastq.gz
    """
}

process downsample_reads {
    
    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory { 1.GB * Math.ceil(Math.max(R1.size(), R2.size()) / 1024 ** 3) }
    time 1.h

    input:
    tuple val(ID), val(LANE), file(R1), file(R2), path('qc-metrics/*', stageAs: './qc-metrics/')
    
    output:
    tuple val(ID), val(LANE), file("${R1.simpleName}.fastq.gz"), file("${R2.simpleName}.fastq.gz"), path("qc-metrics/*")

    script:
    """
    seqkit ID --two-pass -n ${params.read_target} -o ${R1.simpleName}_downsampled.fastq.gz ${R1}
    seqkit ID --two-pass -n ${params.read_target} -o ${R2.simpleName}_downsampled.fastq.gz ${R2}
    mv ${R1.simpleName}_downsampled.fastq.gz ${R1.simpleName}.fastq.gz
    mv ${R2.simpleName}_downsampled.fastq.gz ${R2.simpleName}.fastq.gz
    seqkit stats -j ${task.cpus} -To qc-metrics/downsampling.tsv *.fastq.gz
    """
}