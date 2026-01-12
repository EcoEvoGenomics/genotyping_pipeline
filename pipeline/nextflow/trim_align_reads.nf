#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Trim and align reads
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

include { downsample_reads; deduplicate_reads } from "./qc_utils.nf"
include { qc_reads as qc_raw_reads } from "./qc_utils.nf"
include { qc_reads as qc_trimmed_reads } from "./qc_utils.nf"
include { qc_alignment as qc_unfiltered_alignment } from "./qc_utils.nf"
include { qc_alignment as qc_filtered_alignment } from "./qc_utils.nf"

workflow {
    
    // Take input .CSV with columns ID, LANE, F_READ_PATH, R_READ_PATH
    Channel.fromPath(params.samples)
        .splitCsv()
        .multiMap { cols -> input_reads: [cols[0], cols[1], cols[2], cols[3]] }
        .set { samples }
    
    def input_reads = samples.input_reads
    qc_raw_reads(input_reads, "raw")

    if (params.deduplicate == 'yes') {
        input_reads = deduplicate_reads(input_reads)
    }
    if (params.downsample == 'yes') {
        input_reads = downsample_reads(input_reads)
    }

    def trimmed_reads = trim_reads(input_reads)
    qc_trimmed_reads(trimmed_reads.keys_and_reads, "trimmed")

    // Group separate file pairs (lanes) for each sample and make read groups
    def grouped_reads = trimmed_reads.reads_only \
    | flatten
    | map { file -> 
        def id   = file.name.toString().tokenize("_").get(0)
        return tuple(id, file)
    } \
    | groupTuple(by: 0, sort: true, remainder: true) \
    | group_reads
    
    def unfiltered_alignments = align_reads(grouped_reads, file(params.ref_genome), file(params.ref_index))
    qc_unfiltered_alignment(unfiltered_alignments, file(params.ref_genome), file(params.ref_index), params.ref_scaffold_name, "unfiltered")

    def filtered_alignments = filter_alignment_flags(unfiltered_alignments, file(params.ref_genome), file(params.ref_index), params.exclude_flags)
    qc_filtered_alignment(filtered_alignments, file(params.ref_genome), file(params.ref_index), params.ref_scaffold_name, "filtered")

}

// Step 1 - Read trimming
process trim_reads {

    publishDir "${params.publish_dir}/${ID}/${LANE}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    cpus 4
    memory { task.attempt > 1
        // Double if insufficient in previous attempt - otherwise aim for real peak usage + 50 % from previous
        ? (task.exitStatus == 137 ? task.previousTrace.memory * 2 : task.previousTrace.peak_rss * 1.5)
        // Initial guess
        : 1.MB * Math.max(2048 , 512 * Math.ceil((R1.size() + R2.size()) / 1024 ** 3))
    }
    time { task.attempt > 1 
        // Double if insufficient in previous attempt, otherwise keep previous allocation
        ? (task.exitStatus == 140 ? task.previousTrace.time * 1.5 : task.previousTrace.time )
        // Initial guess
        : 1.m * Math.max(10 , 1.5 * Math.ceil((R1.size() + R2.size()) / 1024 ** 3))
    }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2)

    output:
    tuple val(ID), val(LANE), path("${ID}_${LANE}_R1.fastq.gz"), path("${ID}_${LANE}_R2.fastq.gz"), emit: keys_and_reads
    tuple path("${ID}_${LANE}_R1.fastq.gz"), path("${ID}_${LANE}_R2.fastq.gz"), emit: reads_only

    script:
    """
    fastp \
    --in1 ${R1} \
    --in2 ${R2} \
    --out1 ${ID}_${LANE}_R1.fastq.gz \
    --out2 ${ID}_${LANE}_R2.fastq.gz
    """
}

// Step 3.1 Group read files belonging to same individual (L001, L002, etc.)
process group_reads {
    
    cpus 1
    memory { 8.MB * task.attempt }
    time { 1.m * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), path(grouped_reads, stageAs: "reads/*")

    output:
    tuple val(ID), path(grouped_reads), path("${ID}_reads.list"), env("grouped_reads_total_size")

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
    memory {
        task.attempt == task.maxRetries + 1
        ? 256.GB
        : task.attempt > 1
            ? (task.exitStatus != 140 ? task.previousTrace.memory * 2 : task.previousTrace.memory)
            : 1.GB * Math.max(64, (6 * Math.ceil(grouped_reads_input_size as Long / 1024 ** 3)))
    }
    time {
        task.attempt == task.maxRetries + 1
        ? 12.h
        : task.attempt > 1 
            ? (task.exitStatus == 140 ? task.previousTrace.time * 2 : task.previousTrace.time)
            : 1.m * Math.max(45, (6 * Math.ceil(grouped_reads_input_size as Long / 1024 ** 3)))
    }

    errorStrategy "retry"
    maxRetries 2

    label "require_gpu"

    input:
    tuple val(ID), path(grouped_reads, stageAs: "reads/*"), path(reads_list), val(grouped_reads_input_size)
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
    --out-qc-metrics-dir qc-metrics \
    --tmp .

    # Parse QC files
    rm qc-metrics/*.pdf
    rm qc-metrics/*.png
    for file in qc-metrics/*.txt; do
        mv \$file qc-metrics/${ID}.\$(basename \$file .txt)
    done
    """
}

// Step 4 - Filter alignments by excluding select SAM flags
process filter_alignment_flags {

    publishDir "${params.publish_dir}/${ID}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpus 1
    memory { 1.MB * Math.max(512, 128 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }
    time { 1.m * Math.max(120, 6 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), path(cram), path(index), path(qcmetrics, stageAs: "qc-metrics/*")
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
