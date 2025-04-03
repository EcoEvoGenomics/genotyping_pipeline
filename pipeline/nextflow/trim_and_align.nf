#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Trim and align
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default input parameters
params.ref = file('/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa')
params.publish_dir = './output'

// Workflow
workflow {
    
    // Workflow input
    Channel.fromPath(params.samples)
        .splitCsv()
        .multiMap { cols -> sample: [cols[0], cols[1], cols[2]] }
        .set { samples }

    // Workflow processes
    trim(samples.sample) \
    | align \
    | mark_dup \
    | cram_convert \
    | calc_stats

}

// Step 1 - Read trim
process trim {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input: 
    tuple val(sample), path(f_read), path(r_read)

    output:
    val(sample)
    file("${sample}_R1_TRIM.fastq.gz")
    file("${sample}_R2_TRIM.fastq.gz")
    file("${sample}.html")
    file("${sample}.json")

    script:
    """
    fastp \
    --in1 $f_read \
    --in2 $r_read \
    --out1 ${sample}_R1_TRIM.fastq.gz \
    --out2 ${sample}_R2_TRIM.fastq.gz \
    --report_title "${sample}" \
    --html ${sample}.html \
    --json ${sample}.json
    """
}

// Step 2 - Align to reference genome
process align {

    input:
    val(sample)
    file("${sample}_R1_TRIM.fastq.gz")
    file("${sample}_R2_TRIM.fastq.gz")
    file("${sample}_fastp.html")
    file("${sample}.json")

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

// Step 4 - Mark duplicates in BAM file
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

// Step 5 - Convert BAM file to CRAM file
process cram_convert {
    
    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

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

// Step 6 - Calculate alignment statistics
process calc_stats {

    publishDir "${params.publish_dir}/${sample}/", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    val(sample)
    file(cram)
    file(index)
    file(dedupstats)

    output:
    file("${sample}.cov")
    file("${sample}.flagstat")

    script:
    """
    samtools coverage ${cram} > ${sample}.cov
    samtools flagstat ${cram} > ${sample}.flagstat
    """
}
