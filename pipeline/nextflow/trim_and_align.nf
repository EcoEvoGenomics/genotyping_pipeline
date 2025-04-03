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
    | flatten \
    | map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
    } \
    | groupTuple(by: 0, sort: true, remainder: true) \
    | mark_dup \
    | cram_convert \
    | calc_stats

}

// Step 1 - Read trim
process trim {

    errorStrategy 'ignore'
    publishDir "${params.publish_dir}/trim", saveAs: { filename -> "$filename" }, mode: 'copy'

    input: 
    tuple val(sample), path(f_read), path(r_read)

    output:
    tuple \
    val(sample), \
    path("${sample}_R1_TRIM.fastq.gz"), \
    path("${sample}_R2_TRIM.fastq.gz"), \
    path("${sample}_fastp.html")

    script:
    """
    fastp \
    --in1 $f_read \
    --in2 $r_read \
    --out1 ${sample}_R1_TRIM.fastq.gz \
    --out2 ${sample}_R2_TRIM.fastq.gz \
    --html ${sample}_fastp.html \
    --report-title "${sample}"
    """
}

// Step 2 - Align to reference genome
process align {

    errorStrategy 'ignore'

    input:
    tuple \
    val(sample), \
    path("${sample}_R1_TRIM.fastq.gz"), \
    path("${sample}_R2_TRIM.fastq.gz"), \
    path("${sample}_fastp.html")

    output:
    tuple val(sample), path("${sample}.bam")

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
    echo "Using readgroup: \$READGROUP"

    echo "Aligning ${sample} reads"
    bwa mem -M -t ${task.cpus} -R "\${READGROUP}" ${params.ref} ${sample}_R1_TRIM.fastq.gz ${sample}_R2_TRIM.fastq.gz \
    | samtools view -b | samtools sort -T ${sample} > ${sample}.bam
    """
}

// Step 4 - Mark duplicates in BAM file
process mark_dup {

    errorStrategy 'ignore'

    input:
    tuple val(sample), path("${sample}.bam")

    output:
    tuple val(sample), path("${sample}_dedup.bam"), path("${sample}_dedup.bai")

    script:
    """
    # mark duplicates
    echo "**** Running Picard MarkDuplicates on ${sample} ****"
    picard -Xmx16G MarkDuplicates -I ${sample}.bam -O ${sample}_dedup.bam -M ${sample}_dedup_metrics.txt --TMP_DIR ./run_tmp

    # index bams
    echo "**** Running Picard BuildBamIndex on ${sample} ****"
    picard BuildBamIndex -I ${sample}_dedup.bam --TMP_DIR ./run_tmp
    """
}

// Step 5 - Convert BAM file to CRAM file
process cram_convert {
    
    publishDir "${params.publish_dir}/align", saveAs: { filename -> "$filename" }, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample), path("${sample}_dedup.bam"), path("${sample}_dedup.bai")

    output:
    tuple val(sample), path("${sample}_dedup.cram"), path("${sample}_dedup.cram.crai")

    script:
    """
    samtools view -@ ${task.cpus} -T ${params.ref} -C -o ${sample}_dedup.cram ${sample}_dedup.bam
    samtools index -@ ${task.cpus} ${sample}_dedup.cram 
    """
}

// Step 6 - Calculate alignment statistics
process calc_stats {

    publishDir "${params.publish_dir}/stats", saveAs: { filename -> "$filename" }, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(cram), path(index)

    output:
    tuple path("${sample}_meancov.txt"), path("${sample}.map.stat.csv"), path("${sample}_flagstat.csv")

    script:
    """
    # work out mean coverage
    STATS=\$(samtools depth ${cram} | awk '{sum += \$3} END {print sum / NR}' )
    echo -e "${sample}\t\${STATS}" > ${sample}_meancov.txt

    # run flagstat
    samtools flagstat -@ ${task.cpus} ${cram} > ${sample}_flagstat.csv
    # extract the columns that are wanted
    awk 'OFS="," {print \$1,\$3}' ${sample}_flagstat.csv > ${sample}.col.csv
    # add a header of the sample name
    echo ${sample},${sample} > ${sample}.head.txt
    cat ${sample}.head.txt ${sample}.col.csv > ${sample}.map.stat.csv
    # remove the unneccessary files
    rm ${sample}.col.csv
    rm ${sample}.head.txt
    """
}
