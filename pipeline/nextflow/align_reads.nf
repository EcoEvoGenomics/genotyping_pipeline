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
        def sample_id = dir.getName()
        def r1_file = file("${dir}/${sample_id}_R1_TRIM.fastq.gz")
        def r2_file = file("${dir}/${sample_id}_R2_TRIM.fastq.gz")
        return tuple(sample_id, r1_file, r2_file)
    }
    .set { trimmed_reads }

    def crams = align(trimmed_reads) \
    | mark_dup \
    | cram_convert

    if (params.downsample_crams == 'yes') {
        crams = crams | cram_downsample
    }

    crams | calc_stats
}

// Step 1 - Align to reference genome
process align {

    input:
    tuple \
    val(sample),
    file("${sample}_R1_TRIM.fastq.gz"),
    file("${sample}_R2_TRIM.fastq.gz")

    output:
    val(sample)
    file("${sample}.bam")

    script:
    """
    ### CREATE READ GROUP: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
    FILE_INFO_STRING=\$(zcat ${sample}_R1_TRIM.fastq.gz | head -1)
    INSTRUMENT=\$(echo \${FILE_INFO_STRING} | awk 'BEGIN {FS = ":"}; { print \$1}' | awk '{sub(/@/,""); print}')
    FLOWCELL=\$(echo \${FILE_INFO_STRING} | awk 'BEGIN {FS = ":"}; {print \$3}')
    FLOWCELL_LANE=\$(echo \${FILE_INFO_STRING} | awk 'BEGIN {FS = ":"}; {print \$4}')
    INDEX=\$(echo \${FILE_INFO_STRING} | awk 'BEGIN {FS = ":"}; {print \$10}')
    PLATFORM=Illumina

    FLOWCELL_ID=\${FLOWCELL}.\${FLOWCELL_LANE}
    PLATFORM_UNIT=\${FLOWCELL}.\${FLOWCELL_LANE}.\${INDEX}
    LIBRARY=${sample}.\${INDEX}

    READGROUP="@RG\\tID:\${FLOWCELL_ID}\\tPL:\${PLATFORM}\\tLB:\${LIBRARY}\\tSM:${sample}\\tPU:\${PLATFORM_UNIT}"

    bwa mem -M -t ${task.cpus} -R "\${READGROUP}" ${params.ref} ${sample}_R1_TRIM.fastq.gz ${sample}_R2_TRIM.fastq.gz \
    | samtools view -@ ${task.cpus} -b \
    | samtools sort -@ ${task.cpus} -T ${sample} \
    > ${sample}.bam
    """
}

// Step 2 - Mark duplicates in BAM file
process mark_dup {

    input:
    val(sample)
    file("${sample}.bam")

    output:
    val(sample)
    file("${sample}_dedup.bam")
    file("${sample}_dedup.bai")
    file("${sample}.dedup")

    script:
    """
    picard -Xmx16G MarkDuplicates -I ${sample}.bam -O ${sample}_dedup.bam -M ${sample}.dedup --TMP_DIR ./run_tmp
    picard BuildBamIndex -I ${sample}_dedup.bam --TMP_DIR ./run_tmp
    """
}

// Step 3 - Convert BAM file to CRAM file
process cram_convert {

    input:
    val(sample)
    file("${sample}_dedup.bam")
    file("${sample}_dedup.bai")
    file("${sample}.dedup")

    output:
    val(sample)
    file("${sample}.cram")
    file("${sample}.cram.crai")
    file("${sample}.dedup")

    script:
    """
    samtools view -@ ${task.cpus} -T ${params.ref} -C -o ${sample}.cram ${sample}_dedup.bam
    samtools index -@ ${task.cpus} ${sample}.cram 
    """
}

// Step 4 (Optional) - Downsample CRAM
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

// Step 5 - Calculate alignment statistics
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
