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

process get_alignment_stats {

    publishDir "${params.publish_dir}/${ID}/qc-metrics", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpus 1
    memory { 128.MB * Math.ceil(cram.size() / 1024 ** 3) * task.attempt }
    time { 6.m * Math.ceil(cram.size() / 1024 ** 3) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), file(cram), file(index), file(qcmetrics)
    path(ref_genome)
    path(ref_index)
    val(ref_scaffold_name)
    val(processing_stage)

    output:
    file("${ID}${processing_stage}.tsv")
    file("${ID}${processing_stage}.cramstats")

    script:
    """
    samtools coverage --reference ${ref_genome} ${ID}.cram | grep -v ${ref_scaffold_name} > ${ID}${processing_stage}.tsv
    samtools stats --ref-seq ${ref_genome} ${ID}.cram > ${ID}${processing_stage}.cramstats
    """
}
