process get_reads_stats {

    publishDir "${params.publish_dir}/${ID}/${LANE}/qc-metrics", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    cpus 2
    memory 4.GB
    time 1.h

    input:
    tuple val(ID), val(LANE), file(R1), file(R2)
    val(processing_stage)

    output:
    file("${R1}${processing_stage}.fastqc.zip")
    file("${R2}${processing_stage}.fastqc.zip")

    script:
    """
    zcat ${R1} | fastqc -o ./ stdin:${R1}.fastq.gz
    zcat ${R2} | fastqc -o ./ stdin:${R2}.fastq.gz
    mv ${R1}_fastqc.zip ${ID}_${LANE}_R1${processing_stage}.fastqc.zip
    mv ${R2}_fastqc.zip ${ID}_${LANE}_R2${processing_stage}.fastqc.zip
    """

}

process deduplicate_reads {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory 4.GB
    time 1.h
    
    input:
    tuple val(ID), val(LANE), file(R1), file(R2)

    output:
    tuple val(ID), val(LANE), file("${R1.simpleName}.fastq.gz"), file("${R2.simpleName}.fastq.gz")

    script:
    """
    seqkit rmdup -j ${task.cpus} --by-name -o ${R1.simpleName}_deduplicated.fastq.gz ${R1}
    seqkit rmdup -j ${task.cpus} --by-name -o ${R2.simpleName}_deduplicated.fastq.gz ${R2}
    mv ${R1.simpleName}_deduplicated.fastq.gz ${R1.simpleName}.fastq.gz
    mv ${R2.simpleName}_deduplicated.fastq.gz ${R2.simpleName}.fastq.gz
    """
}

process downsample_reads {
    
    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory { 1.GB * Math.ceil(Math.max(R1.size(), R2.size()) / 1024 ** 3) }
    time 1.h

    input:
    tuple val(ID), val(LANE), file(R1), file(R2)
    
    output:
    tuple val(ID), val(LANE), file("${R1.simpleName}.fastq.gz"), file("${R2.simpleName}.fastq.gz")

    script:
    """
    seqkit ID --two-pass -n ${params.read_target} -o ${R1.simpleName}_downsampled.fastq.gz ${R1}
    seqkit ID --two-pass -n ${params.read_target} -o ${R2.simpleName}_downsampled.fastq.gz ${R2}
    mv ${R1.simpleName}_downsampled.fastq.gz ${R1.simpleName}.fastq.gz
    mv ${R2.simpleName}_downsampled.fastq.gz ${R2.simpleName}.fastq.gz
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
    file("${ID}_${processing_stage}.tsv")
    file("${ID}_${processing_stage}.cramstats")

    script:
    """
    samtools coverage --reference ${ref_genome} ${ID}.cram | grep -v ${ref_scaffold_name} > ${ID}_${processing_stage}.tsv
    samtools stats --ref-seq ${ref_genome} ${ID}.cram > ${ID}_${processing_stage}.cramstats
    """
}
