process qc_reads {

    publishDir "${params.publish_dir}/${ID}/${LANE}/qc-metrics", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    cpus 2
    memory { 1.MB * Math.max(1024, 128 * Math.ceil((R1.size() + R2.size()) / 1024 ** 3)) * task.attempt }
    time { 1.m * Math.max(20, 3 * Math.ceil((R1.size() + R2.size()) / 1024 ** 3)) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2)
    val(processing_stage)

    output:
    path("${ID}_${LANE}_R1_${processing_stage}_fastqc.zip")
    path("${ID}_${LANE}_R2_${processing_stage}_fastqc.zip")

    script:
    """
    mv ${R1} ${ID}_${LANE}_R1_${processing_stage}.fastq.gz
    mv ${R2} ${ID}_${LANE}_R2_${processing_stage}.fastq.gz
    zcat ${ID}_${LANE}_R1_${processing_stage}.fastq.gz | fastqc -o ./ stdin:${ID}_${LANE}_R1_${processing_stage}.fastq.gz
    zcat ${ID}_${LANE}_R2_${processing_stage}.fastq.gz | fastqc -o ./ stdin:${ID}_${LANE}_R2_${processing_stage}.fastq.gz
    mv ${ID}_${LANE}_R1_${processing_stage}_fastqc.zip ${ID}_${LANE}_R1_${processing_stage}_fastqc.zip
    mv ${ID}_${LANE}_R2_${processing_stage}_fastqc.zip ${ID}_${LANE}_R2_${processing_stage}_fastqc.zip
    """

}

process deduplicate_reads {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory { 256.MB * Math.ceil((R1.size() + R2.size()) / 1024 ** 3) * task.attempt }
    time { 3.m * Math.ceil((R1.size() + R2.size()) / 1024 ** 3) * task.attempt }
    
    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2)

    output:
    tuple val(ID), val(LANE), path("${R1.simpleName}.fastq.gz"), path("${R2.simpleName}.fastq.gz")

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
    memory { 256.MB * Math.ceil((R1.size() + R2.size()) / 1024 ** 3) * task.attempt }
    time { 3.m * Math.ceil((R1.size() + R2.size()) / 1024 ** 3) * task.attempt }

    errorStrategy "retry"
    maxRetries 3
    
    input:
    tuple val(ID), val(LANE), path(R1), path(R2)
    
    output:
    tuple val(ID), val(LANE), path("${R1.simpleName}.fastq.gz"), path("${R2.simpleName}.fastq.gz")

    script:
    """
    seqkit ID --two-pass -n ${params.read_target} -o ${R1.simpleName}_downsampled.fastq.gz ${R1}
    seqkit ID --two-pass -n ${params.read_target} -o ${R2.simpleName}_downsampled.fastq.gz ${R2}
    mv ${R1.simpleName}_downsampled.fastq.gz ${R1.simpleName}.fastq.gz
    mv ${R2.simpleName}_downsampled.fastq.gz ${R2.simpleName}.fastq.gz
    """
}

process qc_alignment {

    publishDir "${params.publish_dir}/${ID}/qc-metrics", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpus 1
    memory { 1.MB * Math.max(1024, 128 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }
    time { 1.m * Math.max(20, 6 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), path(cram), path(index), path(qcmetrics)
    path(ref_genome)
    path(ref_index)
    val(ref_scaffold_name)
    val(processing_stage)

    output:
    path("${ID}_${processing_stage}.tsv")
    path("${ID}_${processing_stage}.cramstats")

    script:
    """
    samtools coverage --reference ${ref_genome} ${cram} | grep -v ${ref_scaffold_name} > ${ID}_${processing_stage}.tsv
    samtools stats --ref-seq ${ref_genome} ${cram} > ${ID}_${processing_stage}.cramstats
    """
}
