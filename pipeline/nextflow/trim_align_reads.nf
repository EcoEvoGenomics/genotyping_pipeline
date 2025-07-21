#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Trim and align reads
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

include { downsample_reads; deduplicate_reads } from "./optional.nf"

workflow {
    
    // Input 'samples' is .CSV with columns for ID, LANE, F_READ_PATH, R_READ_PATH
    Channel.fromPath(params.samples)
        .splitCsv()
        .multiMap { cols -> input_reads: [cols[0], cols[1], cols[2], cols[3]] }
        .set { samples }

    // First parse input reads:
    def parsed_reads = parse_input(samples.input_reads)
    if (params.deduplicate == 'yes') {
        parsed_reads = deduplicate_reads(parsed_reads)
    }
    if (params.downsample == 'yes') {
        parsed_reads = downsample_reads(parsed_reads)
    }

    // Then, trim with fastp ...
    def trimmed_reads = trim_reads(parsed_reads)

    // Now, group separate file pairs (lanes) for each sample and make read groups:
    def grouped_reads = trimmed_reads[0] \
    | flatten
    | map { file -> 
        def id   = file.name.toString().tokenize("_").get(0)
        return tuple(id, file)
    } \
    | groupTuple(by: 0, sort: true, remainder: true) \
    | group_reads
    
    // Finally, pass all files of each sample together to the aligner:
    def aligned_reads = align_reads(grouped_reads, file(params.ref_genome), file(params.ref_index))
    def deduplicated_alignments = remove_marked_duplicates(aligned_reads, file(params.ref_genome), file(params.ref_index), params.exclude_flags)
    get_final_alignment_stats(deduplicated_alignments, file(params.ref_genome), file(params.ref_index), params.ref_scaffold_name)

}

// Step 1 - Parse input reads for a sample
process parse_input {

    container "quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
    cpus 2
    memory { 1.MB * Math.max(128 , 12 * Math.ceil((R1.size() + R2.size())) / 1024 ** 3) * task.attempt }
    time { 1.m * Math.max(10 , 1 * Math.ceil((R1.size() + R2.size())) / 1024 ** 3) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2)

    output:
    tuple val(ID), val(LANE), path(R1), path(R2), path('qc-metrics/*')

    script:
    """
    mkdir qc-metrics
    seqkit stats -j ${task.cpus} -To qc-metrics/unprocessed.tsv *.fastq.gz
    printf '%s\\t%s\\n' \
        'ID_LANE' '${ID}_${LANE}' \
        'deduplicate' '${params.deduplicate}' \
        'downsample' '${params.downsample}' \
        > qc-metrics/settings.tsv
    """
}

// Step 2 - Read trim_reads
process trim_reads {

    publishDir "${params.publish_dir}/${ID}/${LANE}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    cpus 4
    memory { 1.MB * Math.max(1024 , 256 * Math.ceil((R1.size() + R2.size())) / 1024 ** 3) * task.attempt }
    time { 1.m * Math.max(10 , 1 * Math.ceil((R1.size() + R2.size())) / 1024 ** 3) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), file(R1), file(R2), path(qcmetrics, stageAs: './qc-metrics/')

    output:
    tuple file("${ID}_${LANE}_R1.fastq.gz"), file("${ID}_${LANE}_R2.fastq.gz")
    path("qc-metrics/*")

    script:
    """
    fastp \
    --in1 ${R1} \
    --in2 ${R2} \
    --out1 ${ID}_${LANE}_R1.fastq.gz \
    --out2 ${ID}_${LANE}_R2.fastq.gz \
    --report_title "${ID}_${LANE}" \
    --html qc-metrics/${ID}_${LANE}.html \
    --json qc-metrics/${ID}_${LANE}.json
    """
}

// Step 3.1 Group read files belonging to same individual (L001, L002, etc.)
process group_reads {
    
    cpus 1
    memory { 4.MB * task.attempt }
    time { 10.s * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), path(grouped_reads, stageAs: "reads/*")

    output:
    tuple val(ID), path(grouped_reads), file("${ID}_reads.list"), env("grouped_reads_total_size")

    script:
    """
    # Obtain list of lanes present for the sample
    find -type l,f -name "${ID}*.fastq.gz" \
    | awk -F'_' '{print \$2}' \
    | sort \
    | uniq \
    > lanes.list

    # For each lane, write entry to fastq file list
    while IFS= read -r lane; do

        # Obtain info for for read group (https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
        READ_HEADER=\$(zcat reads/${ID}_\${lane}_R1.fastq.gz | head -1)

        INSTRUMENT=\$(echo \${READ_HEADER} | awk 'BEGIN {FS = ":"}; { print \$1}' | awk '{sub(/@/,""); print}')
        FLOWCELL=\$(echo \${READ_HEADER} | awk 'BEGIN {FS = ":"}; {print \$3}')
        INDEX=\$(echo \${READ_HEADER} | awk 'BEGIN {FS = ":"}; {print \$10}')
        PLATFORM=Illumina

        # Construct read group
        FLOWCELL_ID=\${FLOWCELL}.\${lane}
        PLATFORM_UNIT=\${FLOWCELL_ID}.\${INDEX}
        LIBRARY=${ID}.\${INDEX}
        
        READGROUP="@RG\\tID:\${FLOWCELL_ID}\\tLB:\${LIBRARY}\\tPL:\${PLATFORM}\\tSM:${ID}\\tPU:\${PLATFORM_UNIT}"

        # Write to list
        printf '%s %s %s\\n' \
            "reads/${ID}_\${lane}_R1.fastq.gz" \
            "reads/${ID}_\${lane}_R2.fastq.gz" \
            "\${READGROUP}" \
        >> ${ID}_reads.list

    done < lanes.list

    # Calculate total size of grouped read files
    grouped_reads_total_size=\$( du -sbL reads/ | awk '{ print \$1 }' )
    """
}

// Step 3.2 - Align to reference genome using GPU and obtain stats pre-removal of marked duplicates
process align_reads {

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--nv"
    memory { 1.GB * Math.max(36, 8 + (4 * Math.ceil(grouped_reads_input_size as Long / 1024 ** 3))) * (1 + (0.25 * (task.attempt - 1))) }
    time { 1.m * Math.max(15, 10 + (1 * Math.ceil(grouped_reads_input_size as Long / 1024 ** 3))) * (1 + (0.25 * (task.attempt - 1))) }

    errorStrategy "retry"
    maxRetries 3

    label "require_gpu"

    input:
    tuple val(ID), path(grouped_reads, stageAs: "reads/*"), file(reads_list), val(grouped_reads_input_size)
    path(reference_genome)
    path(reference_genome_index)

    output:
    tuple val(ID), path("${ID}_marked.cram"), path("${ID}_marked.cram.crai"), path('qc-metrics/*')

    script:
    """
    # Prepare QC directory
    mkdir qc-metrics

    # Align with Parabricks FQ2BAM
    pbrun fq2bam \
    --ref ${reference_genome} \
    --in-fq-list ${reads_list} \
    --out-bam ${ID}_marked.cram \
    --out-duplicate-metrics qc-metrics/dedup.txt \
    --out-qc-metrics-dir qc-metrics \
    --tmp .

    # Obtain coverage statistics with Parabricks BAMMETRICS
    pbrun bammetrics \
    --ref ${reference_genome} \
    --bam ${ID}_marked.cram \
    --out-metrics-file qc-metrics/bammetrics.txt \
    --tmp-dir .

    # Parse QC files
    rm qc-metrics/*.pdf
    rm qc-metrics/*.png
    for file in qc-metrics/*.txt; do
        mv \$file qc-metrics/${ID}.\$(basename \$file .txt)
    done
    """
}

// Step 4 - Remove duplicates marked in the previous process
process remove_marked_duplicates {

    publishDir "${params.publish_dir}/${ID}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpus 1
    memory { 1.MB * Math.max(1024, 128 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }
    time { 1.m * Math.max(15, 3 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), file(cram), file(index), path(qcmetrics, stageAs: "qc-metrics/*")
    path(ref_genome)
    path(ref_index)
    val(exclude_flags)

    output:
    tuple val(ID), path("${ID}.cram"), path("${ID}.cram.crai"), path('qc-metrics/*', includeInputs: true)

    script:
    """
    samtools view -@ ${task.cpus} --reference ${ref_genome} --cram --excl-flags ${exclude_flags} ${cram} > ${ID}.cram
    samtools index -@ ${task.cpus} ${ID}.cram
    """
}

// Step 4 - Obtain stats for final alignment
process get_final_alignment_stats {

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

    output:
    file("${ID}.tsv")
    file("${ID}.cramstats")

    script:
    """
    samtools coverage --reference ${ref_genome} ${ID}.cram | grep -v ${ref_scaffold_name} > ${ID}.tsv
    samtools stats --ref-seq ${ref_genome} ${ID}.cram > ${ID}.cramstats
    """
}
