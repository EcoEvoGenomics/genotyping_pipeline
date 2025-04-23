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
    memory { 8.GB + 10.GB * Math.ceil(r1_file.size() / 1024 ** 3) }

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
