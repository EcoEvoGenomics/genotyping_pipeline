include { downsample_reads; deduplicate_reads } from "./qc_utils.nf"
include { qc_reads as qc_raw_reads } from "./qc_utils.nf"
include { qc_reads as qc_trimmed_reads } from "./qc_utils.nf"
include { qc_alignment as qc_unfiltered_alignment } from "./qc_utils.nf"
include { qc_alignment as qc_filtered_alignment } from "./qc_utils.nf"

workflow {
    
    // Take input .CSV with columns ID, SEX, LANE, F_READ_PATH, R_READ_PATH
    Channel.fromPath(params.samples)
        .splitCsv()
        .multiMap { cols -> input_reads: [cols[0], cols[2], cols[3], cols[4]] }
        .set { samples }

    // Fetch reference index files
    def ref_index = file(params.ref_genome.toString() + ".fai")
    def ref_indices = files(params.ref_genome.toString() + "*.{amb,ann,bwt,fai,pac,sa}")
    
    // Read preprocessing
    def input_reads = samples.input_reads
    if (params.deduplicate) {
        input_reads = deduplicate_reads(input_reads)
    }
    if (params.downsample) {
        input_reads = downsample_reads(input_reads)
    }
    def trimmed_reads = trim_reads(input_reads, file(params.trimming_flags))
    def readgrouped_trimmed_reads = append_readgroups(trimmed_reads)
    
    // Read quality control
    qc_raw_reads(input_reads, "raw")
    qc_trimmed_reads(trimmed_reads, "trimmed")

    // Alignment
    if (params.aligner != 'gpu') {
        if (params.aligner == "mem") {
            per_lane_unfiltered_alignments = align_bwa_mem(
                readgrouped_trimmed_reads.bwa_input,
                file(params.ref_genome),
                ref_indices
            )
        }
        if(params.aligner == "aln") {
            per_lane_unfiltered_alignments = align_bwa_aln(
                readgrouped_trimmed_reads.bwa_input,
                file(params.ref_genome),
                ref_indices
            )
        }
        grouped_unfiltered_alignments = per_lane_unfiltered_alignments \
        | flatten \
        | map { file -> 
            def sample_id = file.name.toString().tokenize("_").get(0)
            return tuple(sample_id, file)
        } \
        | groupTuple(by: 0, sort: true, remainder: true)
        merged_unfiltered_alignments = merge_lanes(grouped_unfiltered_alignments, params.ref_genome, ref_indices)
        sorted_unfiltered_alignments = sort_alignment(merged_unfiltered_alignments, params.ref_genome, ref_indices)
        unfiltered_alignments = mark_duplicates(sorted_unfiltered_alignments, params.ref_genome, ref_indices)
    } else {
        def reads_for_gpu = readgrouped_trimmed_reads.fq2bam_input \
        | flatten \
        | map { file ->
            def sample_id = file.name.toString().tokenize("_").get(0)
            return tuple(sample_id, file)
        } 
        | groupTuple(by: 0, sort: true, remainder: true) \
        | parse_input_for_fq2bam
        unfiltered_alignments = align_fq2bam(reads_for_gpu, file(params.ref_genome), ref_index)
    }
    def filtered_alignments = filter_alignment(unfiltered_alignments, file(params.ref_genome), ref_index, params.exclude_flags)

    // Alignment quality control
    qc_unfiltered_alignment(unfiltered_alignments, file(params.ref_genome), ref_index, params.ref_scaffold_name, "unfiltered")
    qc_filtered_alignment(filtered_alignments, file(params.ref_genome), ref_index, params.ref_scaffold_name, "filtered")

}

process trim_reads {

    publishDir "${params.publish_dir}/${ID}/${LANE}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/fastp:0.24.0--heae3180_1"
    cpus 4
    memory { task.attempt > 1
        // Double if insufficient in previous attempt - otherwise aim for real peak usage + 50 % from previous
        ? (task.exitStatus == 137 ? task.previousTrace.memory * 2 : task.previousTrace.peak_rss * 1.5)
        // Initial guess
        : 1.MB * Math.max(8192 , 512 * Math.ceil((R1.size() + R2.size()) / 1024 ** 3))
    }
    time { task.attempt > 1 
        // Double if insufficient in previous attempt, otherwise keep previous allocation
        ? (task.exitStatus == 140 ? task.previousTrace.time * 2 : task.previousTrace.time )
        // Initial guess
        : 3.m * Math.max(30 , 1.5 * Math.ceil((R1.size() + R2.size()) / 1024 ** 3))
    }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2)
    path(filterfile)

    output:
    tuple val(ID), val(LANE), path("${ID}_${LANE}_R1.fastq.gz"), path("${ID}_${LANE}_R2.fastq.gz")

    script:
    """
    echo '--in1 ${R1}' >> fastp.args
    echo '--in2 ${R2}' >> fastp.args
    echo '--out1 ${ID}_${LANE}_R1.fastq.gz' >> fastp.args
    echo '--out2 ${ID}_${LANE}_R2.fastq.gz' >> fastp.args
    cat ${filterfile} >> fastp.args

    cat fastp.args | xargs fastp
    """
}

process append_readgroups {

    cpus 1
    memory { 8.MB * task.attempt }
    time { 1.m * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2)

    output:
    tuple val(ID), val(LANE), path("${ID}_${LANE}_R1.fastq.gz"), path("${ID}_${LANE}_R2.fastq.gz"), env("READGROUP"), emit: bwa_input
    tuple path("${ID}_${LANE}_R1.fastq.gz"), path("${ID}_${LANE}_R2.fastq.gz"), path("${ID}_${LANE}_GROUPED.txt"), emit: fq2bam_input

    script:
    """
    # Obtain info for for read group (https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
    READ_HEADER=\$(zcat ${ID}_${LANE}_R1.fastq.gz | head -1)

    INSTRUMENT=\$(echo \${READ_HEADER} | awk 'BEGIN {FS = ":"}; { print \$1}' | awk '{sub(/@/,""); print}')
    FLOWCELL=\$(echo \${READ_HEADER} | awk 'BEGIN {FS = ":"}; {print \$3}')
    INDEX=\$(echo \${READ_HEADER} | awk 'BEGIN {FS = ":"}; {print \$10}')
    PLATFORM=Illumina

    # Construct read group
    FLOWCELL_ID=\${FLOWCELL}.${LANE}
    PLATFORM_UNIT=\${FLOWCELL_ID}.\${INDEX}
    LIBRARY=${ID}.\${INDEX}
    
    READGROUP="@RG\\tID:\${FLOWCELL_ID}\\tLB:\${LIBRARY}\\tPL:\${PLATFORM}\\tSM:${ID}\\tPU:\${PLATFORM_UNIT}"
    printf '%s %s %s\\n' \
        "reads/${ID}_${LANE}_R1.fastq.gz" \
        "reads/${ID}_${LANE}_R2.fastq.gz" \
        "\${READGROUP}" \
    >> ${ID}_${LANE}_GROUPED.txt
    """
}

process parse_input_for_fq2bam {
    
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
    find . -type l,f -name "${ID}*_GROUPED.txt" > lane_readgroup_files.txt
    while IFS= read -r lane_readgroup_file; do
        cat \$lane_readgroup_file >> ${ID}_reads.list
    done < lane_readgroup_files.txt
    grouped_reads_total_size=\$( du -sbL reads/ | awk '{ print \$1 }' )
    """
}

process align_fq2bam {

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
    path(ref_genome)
    path(ref_index)

    output:
    tuple val(ID), path("${ID}_marked.cram"), path("${ID}_marked.cram.crai"), path('qc-metrics/*')

    script:
    """
    # Prepare QC directory
    mkdir qc-metrics

    # Align with Parabricks FQ2BAM
    pbrun fq2bam \
    --ref ${ref_genome} \
    --in-fq-list ${reads_list} \
    --bwa-options="-K 10000000" \
    --out-bam ${ID}_marked.cram \
    --out-qc-metrics-dir qc-metrics \
    --out-duplicate-metrics qc-metrics/duplicates.txt \
    --tmp .

    # Parse QC files
    rm qc-metrics/*.pdf
    rm qc-metrics/*.png
    for file in qc-metrics/*.txt; do
        mv \$file qc-metrics/${ID}.\$(basename \$file .txt)
    done
    """
}

process align_bwa_mem {
    
    // Container build page: https://wave.seqera.io/view/builds/bd-7673b9c618b4118a_1
    container "community.wave.seqera.io/library/bwa_samtools:7673b9c618b4118a"
    cpus 16
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2), val(readgroup)
    path(ref_genome)
    path(ref_and_bwa_indices, stageAs: "./*")

    output:
    path("${ID}_${LANE}.cram")

    script:
    """
    bwa mem -t ${task.cpus} -K 10000000 -R '${readgroup}' ${ref_genome} ${R1} ${R2} \
    | samtools view -@ ${task.cpus} -T ${ref_genome} -C | samtools sort -T ${ID} -o ${ID}_${LANE}.cram
    """
}

process align_bwa_aln {
    
    // Container build page: https://wave.seqera.io/view/builds/bd-7673b9c618b4118a_1
    container "community.wave.seqera.io/library/bwa_samtools:7673b9c618b4118a"
    cpus 16
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), val(LANE), path(R1), path(R2), val(readgroup)
    path(ref_genome)
    path(ref_and_bwa_indices, stageAs: "./*")

    output:
    path("${ID}_${LANE}.cram")

    script:
    """
    bwa aln -t ${task.cpus} ${ref_genome} ${R1} > R1.sai
    bwa aln -t ${task.cpus} ${ref_genome} ${R2} > R2.sai
    bwa sampe -r '${readgroup}' ${ref_genome} R1.sai R2.sai ${R1} ${R2} \
    | samtools view -@ ${task.cpus} -T ${ref_genome} -C | samtools sort -T ${ID} -o ${ID}_${LANE}.cram
    """
}

process merge_lanes {
    
    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(ID), path(crams, stageAs: "crams/*")
    path(ref_genome)
    path(ref_and_bwa_indices, stageAs: "./*")

    output:
    tuple val(ID), path("${ID}.cram")

    script:
    def cramlist = crams instanceof List ? crams.join(" ") : crams
    """
    samtools merge -@ ${task.cpus} -rf ${ID}.cram ${cramlist}
    """
}

process sort_alignment {
    
    container "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
    cpus 1
    memory { 32.GB * task.attempt }
    time { 8.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    label "high_mem_per_cpu"

    input:
    tuple val(ID), path(cram)
    path(ref_genome)
    path(ref_and_bwa_indices, stageAs: "./*")

    output:
    tuple val(ID), path("${ID}_sorted.cram")

    script:
    """
    gatk SortSam \
        --java-options -Xmx${task.memory.toGiga()}G \
        --MAX_RECORDS_IN_RAM 5000000 \
        --TMP_DIR . \
        -I ${cram} \
        -O ${ID}_sorted.cram \
        -R ${ref_genome} \
        --SORT_ORDER coordinate
    """
}

process mark_duplicates {

    // Container build page: https://wave.seqera.io/view/builds/bd-77c6fcf88cba7ceb_1
    container "community.wave.seqera.io/library/gatk4_samtools:77c6fcf88cba7ceb"
    cpus 1
    memory { 32.GB * task.attempt }
    time { 8.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    label "high_mem_per_cpu"

    input:
    tuple val(ID), path(cram)
    path(ref_genome)
    path(ref_and_bwa_indices, stageAs: "./*")

    output:
    tuple val(ID), path("${ID}_marked.cram"), path("${ID}_marked.cram.crai"), path('qc-metrics/*')

    script:
    """
    mkdir qc-metrics
    
    gatk MarkDuplicates \
        --java-options -Xmx${task.memory.toGiga()}G \
        --TMP_DIR . \
        -I ${cram} \
        -O ${ID}_marked.cram \
        -R ${ref_genome} \
        -M qc-metrics/duplicates.txt
    samtools index -@ ${task.cpus} ${ID}_marked.cram
    
    for file in qc-metrics/*.txt; do
        mv \$file qc-metrics/${ID}.\$(basename \$file .txt)
    done
    """
}

process filter_alignment {

    publishDir "${params.publish_dir}/${ID}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpus 1
    memory { 1.MB * Math.max(2048, 128 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }
    time { 1.m * Math.max(240, 6 * Math.ceil(cram.size() / 1024 ** 3)) * task.attempt }

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
