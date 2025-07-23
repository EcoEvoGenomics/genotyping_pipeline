#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Phase VCF
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

include { define_windows } from "./call_variants.nf"

workflow{

    define_windows(params.ref_index, params.window_size, params.ref_scaffold_name)
    def window_list = define_windows.out.windows.map{path -> file(path.toString())}.readLines()
    def unphased_vcf = Channel.fromPath("${params.unphased_vcf}")
    def unphased_csi = Channel.fromPath("${params.unphased_csi}")
    
    phase_common_variants(window_list, params.ref_recombination_map_dir, params.ref_scaffold_name, unphased_vcf, unphased_csi) \
    | flatten \
    | map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
    } \
    | groupTuple( by:0, sort:true ) \
    | ligate_phased \
    | convert_bcf_to_vcf

}

process phase_common_variants {
    
    container "quay.io/biocontainers/shapeit5:5.1.1--hb60d31d_0"
    cpus 4
    memory { 1.GB + (1.GB * Math.log(Math.ceil(unphased_vcf.size() / 1024 ** 3))) * task.attempt }
    time { 1.m * Math.ceil(unphased_vcf.size() / 1024 ** 3) * task.attempt }

    // SHAPEIT5 exits with an error if there are no variants in the window. If so, ignore the window.
    errorStrategy { task.exitStatus in 137..140 ? "retry" : "ignore" }
    maxRetries 3

    input:
    each(window)
    val(recombination_map_dir)
    val(ref_scaffold_name)
    path(unphased_vcf)
    path(unphased_csi)

    output:
    tuple path("${window}.bcf"), path("${window}.bcf.csi"), optional: true

    script:
    """
    WINDOW=${window}
    
    if [[ \$WINDOW == *"${ref_scaffold_name}"* ]]; then
        exit 0;
    fi
    
    MAP="${recombination_map_dir}/\${WINDOW%:*}.map"
    if [[ -f \${MAP} ]]; then
        SHAPEIT5_phase_common --input ${unphased_vcf} --region ${window} --map \${MAP} --output ${window}.bcf --thread ${task.cpus}
    else
        SHAPEIT5_phase_common --input ${unphased_vcf} --region ${window} --output ${window}.bcf --thread ${task.cpus}
    fi
    """
}

process ligate_phased {

    publishDir "${params.publish_dir}/${key}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/shapeit5:5.1.1--hb60d31d_0"
    cpus 4
    memory { 1.GB * task.attempt }
    time { 30.m * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(key), path(bcfs_with_csis)

    output:
    tuple val(key), path("${key}.bcf"), path("${key}.bcf.csi") 

    script:
    """
    sorted_bcfs=\$(echo ${bcfs_with_csis} | tr ' ' '\n' | grep -v ".csi" | sort -t"-" -k2 -n | tr '\n' ' ')
    for i in \$sorted_bcfs; do echo \$i >> chunks.txt ; done
    SHAPEIT5_ligate --input chunks.txt --output ${key}.bcf --index --thread ${task.cpus}
    """
}

process convert_bcf_to_vcf {

    publishDir "${params.publish_dir}/${key}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 2
    memory 1.GB
    time 2.h

    input:
    tuple val(key), path(bcf), path(csi)

    output:
    tuple val(key), path("${key}.vcf.gz"), path("${key}.vcf.gz.csi")

    script:
    """
    bcftools view --output-type z --output ${key}.vcf.gz ${bcf}
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """

}
