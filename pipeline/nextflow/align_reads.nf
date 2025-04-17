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

    align_gpu(params.ref, params.ref_index, trimmed_reads)

    // if (params.downsample_crams == 'yes') {
    //     crams = crams | cram_downsample
    // }

    // crams | calc_stats
}

// Step 1 - Align to reference genome using GPU
process align_gpu {
    
    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    path('ref_genome.fa')
    path('ref_genome.fa.fai')
    tuple \
    val(sample), \
    path(fwd, stageAs: "R1.fastq.gz"), \
    path(rev, stageAs: "R2.fastq.gz")

    output:
    val(sample)
    path("${sample}.cram")
    path("${sample}.cram.crai")
    path("${sample}.dedup")
    path('qc-metrics/*')

    script:
    """
    pbrun fq2bam \
    --ref ref_genome.fa \
    --in-fq R1.fastq.gz R2.fastq.gz \
    --out-bam ${sample}.cram \
    --out-duplicate-metrics ${sample}.dedup \
    --out-qc-metrics-dir qc-metrics \
    --tmp .
    """
}

// Step 2 (Optional) - Downsample CRAM
process cram_downsample {
    
    input:
    val(sample)
    file('original.cram')
    file('original.cram.crai')
    file(dedupstats)

    output:
    val(sample)
    file("${sample}.cram")
    file("${sample}.cram.crai")
    file("${sample}.dedup")

    script:
    """
    mean_coverage=\$(samtools depth original.cram | awk '{sum += \$3} END {result = sum / NR; printf "%d\\n", result}' )

    if (( \$(echo "scale=4; \${mean_coverage} > ${params.max_cram_depth}" | bc -l) )); then
        fraction_sampled=\$(echo "scale=4; ${params.max_cram_depth} / \${mean_coverage}" | bc)
        samtools view -@ ${task.cpus} -s \${fraction_sampled} -C -T ${params.ref} -o ${sample}.cram original.cram
        samtools index -@ ${task.cpus} ${sample}.cram
    else
        mv original.cram ${sample}.cram
        mv original.cram.crai ${sample}.cram.crai
    fi
    """
}

// Step 3 - Calculate alignment statistics
process calc_stats {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    val(sample)
    file(cram)
    file(index)
    file(dedupstats)

    output:
    file("${sample}.dedup")
    file("${sample}.cram")
    file("${sample}.cram.crai")
    file("${sample}.tsv")
    file("${sample}.flagstat")

    script:
    """
    samtools coverage ${cram} | grep -v ${params.ref_scaffold_name} > ${sample}.tsv
    samtools flagstat ${cram} > ${sample}.flagstat
    """
}
