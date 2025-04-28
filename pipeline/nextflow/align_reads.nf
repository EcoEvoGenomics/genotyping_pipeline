#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Align reads
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

workflow {
    
    Channel.fromPath("${params.reads_dir}/*", type: 'dir')
    .map { dir -> 
        def sample = dir.getName()
        def r1_file = file("${dir}/${sample}_R1_TRIM.fastq.gz")
        def r2_file = file("${dir}/${sample}_R2_TRIM.fastq.gz")
        return tuple(sample, r1_file, r2_file)
    }
    .set { trimmed_reads }

    def aligned = align_reads(trimmed_reads, params.ref_genome, params.ref_index)
    calc_coverage(aligned, params.ref_genome, params.ref_index)
}

// Step 1 - Align to reference genome using GPU
process align_reads {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--nv"
    
    clusterOptions "--job-name=align_reads --gpus=2"
    cpus 24
    memory { 16.GB + 7.5.GB * Math.ceil(Math.max(r1_file.size(), r2_file.size()) / 1024 ** 3) }
    time 1.h
    
    label "gpu"

    input:
    tuple val(sample), path(r1_file), path(r2_file)
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

// Step 2 - Additional alignment statistics
process calc_coverage {

    publishDir "${params.publish_dir}/${sample}/qc-metrics", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"

    clusterOptions "--job-name=calc_coverage"
    cpus 1
    memory 1.GB
    time 1.h

    input:
    tuple \
    val(sample), \
    file(cram), \
    file(index), \
    file(qcmetrics)
    path(ref_genome)
    path(ref_index)

    output:
    file("${sample}.tsv")
    file("${sample}.flagstat")

    script:
    """
    samtools coverage --reference ${ref_genome} ${cram} | grep -v ${params.ref_scaffold_name} > ${sample}.tsv
    samtools flagstat ${cram} > ${sample}.flagstat
    """
}
