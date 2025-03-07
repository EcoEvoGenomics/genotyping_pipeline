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
params.trim = file('/cluster/projects/nn10082k/trimmomatic_adapters/')
params.publish_dir = './output'

// Workflow
workflow{
    Channel
    .from(file(params.samples))
    .splitCsv()
    .multiMap { it ->
        samples: [it[0]]
        new_samples: [it[1]]
        f_reads: [it[0], it[1], it[2]]
        r_reads: [it[0], it[1], it[3]]
        adapter: it[4]
        }
    .set{result}

    trimming(result.f_reads, result.r_reads, result.adapter) | align \
    // this will take all the outputs and group them - accounting for other lanes
    | flatten \
    | map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
    } \
    | groupTuple(by: 0, sort: true, remainder: true) \
    | merge_sort \
    | mark_dup | cram_convert | calc_stats 
}

// Step 1 - Read trimming
process trimming {

    errorStrategy 'ignore'
    publishDir '${params.publish_dir}/trim', saveAs: { filename -> "$filename" }

    input: 
    tuple val(f_sample), val(f_new_sample), path(f_read)
    tuple val(sample), val(new_sample), path(r_read)
    val(adapter)

    output:
    tuple \
    val(new_sample), \
    path(f_read), \
    path(r_read), \
    path("${new_sample}.R1.trim_pair.fastq.gz"), \
    path("${new_sample}.R2.trim_pair.fastq.gz"), \
    path("${new_sample}.R1.trim_unpair.fastq.gz"), \
    path("${new_sample}.R2.trim_unpair.fastq.gz"), \
    path("${new_sample}.stats"), \
    path("${f_read}_fastqc.zip"), \
    path("${r_read}_fastqc.zip")

    script:
    """
    ## run fastqc
    zcat ${f_read} | fastqc -o ./ stdin:${f_read}.fastq.gz
    zcat ${r_read} | fastqc -o ./ stdin:${r_read}.fastq.gz

    ## set the adapter fasta - need to find a way to change this
    ADAPT_FAST=${params.trim}/${adapter}.fa
    ## run trimmometic
    trimmomatic PE -threads ${task.cpus} $f_read $r_read \
    ${new_sample}.R1.trim_pair.fastq.gz ${new_sample}.R1.trim_unpair.fastq.gz \
    ${new_sample}.R2.trim_pair.fastq.gz ${new_sample}.R2.trim_unpair.fastq.gz \
    ILLUMINACLIP:\${ADAPT_FAST}:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:10 MINLEN:50 \
    |& tee ${new_sample}.stats
    """
}

// Step 2 - Align to reference genome
process align {

    errorStrategy 'ignore'

    input:
    tuple \
    val(sample), \
    path(f_read), \
    path(r_read), \
    path("${sample}.R1.trim_pair.fastq.gz"), \
    path("${sample}.R2.trim_pair.fastq.gz"), \
    path("${sample}.R1.trim_unpair.fastq.gz"), \
    path("${sample}.R2.trim_unpair.fastq.gz"), \
    path("${sample}.stats"), \
    path("${f_read}_fastqc.zip"), \
    path("${r_read}_fastqc.zip")

    output:
    tuple \
    path("${sample}_pair.bam"), \
    path("${sample}_F_unpair.bam"), \
    path("${sample}_R_unpair.bam")

    script:
    """
    ### CREATE READ GROUP
    # create base string from file info
    STRING=\$(zcat ${sample}.R1.trim_pair.fastq.gz | head -1)
    # break string into information
    # instrument
    INSTRUMENT=\$(echo \${STRING} | awk 'BEGIN {FS = ":"}; { print \$1}' | awk '{sub(/@/,""); print}')
    # run id
    RUN=\$(echo \${STRING} | awk 'BEGIN {FS = ":"}; {print \$2}')
    # flowcell
    FLOWCELL=\$(echo \${STRING} | awk 'BEGIN {FS = ":"}; {print \$3}')
    # flowcell lane
    FLOWCELL_LANE=\$(echo \${STRING} | awk 'BEGIN {FS = ":"}; {print \$4}')
    # index sequence
    INDEX=\$(echo \${STRING} | awk 'BEGIN {FS = ":"}; {print \$10}')
    # platform (always Illumina)
    PLATFORM=Illumina

    ## construct read group - see here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
    # first construct ID - must be unique - flowcell and number
    ID=\${FLOWCELL}.\${FLOWCELL_LANE}
    # next PU (platform unit) - flowcell , lane and sample barcode
    PU_DATA=\${FLOWCELL}.\${FLOWCELL_LANE}.\${INDEX}
    # set library
    LIBRARY=${sample}.\${INDEX}
    # final read group config
    READGROUP="@RG\\tID:\${ID}\\tPL:\${PLATFORM}\\tLB:\${LIBRARY}\\tSM:${sample}\\tPU:\${PU_DATA}"

    echo "Using readgroup: \$READGROUP"

    ### MAP PAIRED
    echo "Aligning ${sample} paired reads"
    # run the alignment on paired reads
    bwa mem -M -t ${task.cpus} -R "\${READGROUP}" ${params.ref} ${sample}.R1.trim_pair.fastq.gz ${sample}.R2.trim_pair.fastq.gz \
    | samtools view -b | samtools sort -T ${sample} > ${sample}_pair.bam


    ### MAP UNPAIR FORWARD
    echo "Aligning ${sample} unpaired forward reads."
    # run alignment
    bwa mem -M -t ${task.cpus} -R "\${READGROUP}" ${params.ref} ${sample}.R1.trim_unpair.fastq.gz \
    | samtools view -b | samtools sort -T ${sample} > ${sample}_F_unpair.bam


    ### MAP UNPAIR REVERSE
    echo "Aligning ${sample} unpaired reverse reads."
    # run alignment
    bwa mem -M -t ${task.cpus} -R "\${READGROUP}" ${params.ref} ${sample}.R2.trim_unpair.fastq.gz \
    | samtools view -b | samtools sort -T ${sample} > ${sample}_R_unpair.bam

    """
}

// Step 3 - Merge and sort BAM file
process merge_sort {

    errorStrategy 'ignore'

    input:
    tuple val(sample), path(bams, stageAs: "?/bam?.bam")

    output:
    tuple val(sample), path("${sample}_merge_sort.bam")

    
    script:
    def bam_list = bams instanceof List ? bams.join(" ") : bams
    """
    ### MERGE SAMPLES
    echo "Merging bams for ${sample}"
    # merge
    samtools merge -rf -@ ${task.cpus} ${sample}_merge.bam ${bam_list}
    # sort
    echo "Sorting merged bam for ${sample}"
    samtools sort -@ ${task.cpus} -T ${sample}_tmp -o ${sample}_merge_sort.bam ${sample}_merge.bam
    """
}

// Step 4 - Mark duplicates in BAM file
process mark_dup {

    errorStrategy 'ignore'

    input:
    tuple val(sample), path("merge_sort.bam")

    output:
    tuple val(sample), path("${sample}_dedup.bam"), path("${sample}_dedup.bai")

    script:
    """
    # mark duplicates
    echo "**** Running Picard MarkDuplicates on ${sample} ****"
    picard -Xmx16G MarkDuplicates -I merge_sort.bam -O ${sample}_dedup.bam -M ${sample}_dedup_metrics.txt --TMP_DIR ./run_tmp

    # index bams
    echo "**** Running Picard BuildBamIndex on ${sample} ****"
    picard BuildBamIndex -I ${sample}_dedup.bam --TMP_DIR ./run_tmp
    """
}

// Step 5 - Convert BAM file to CRAM file
process cram_convert {
    
    publishDir '${params.publish_dir}/align', saveAs: { filename -> "$filename" }
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

    publishDir '${params.publish_dir}/stats', saveAs: { filename -> "$filename" }, mode: 'copy'
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
