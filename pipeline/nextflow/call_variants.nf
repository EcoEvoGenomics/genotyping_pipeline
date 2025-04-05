#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Call VCF
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Workflow
workflow{    

    windows_list = file("${params.windows_dir}/genome_windows.list").readLines()
        
    Channel.fromPath(params.ploidy_file).set{ploidy_file}
    Channel.fromPath(params.windows_dir).set{windows_dir}
    Channel.fromList(windows_list).set{windows}
    Channel.from(file(params.crams)).set{crams}

    def chromosome_vcfs = genotype_windows(crams, ploidy_file, windows_dir, windows) \
    | map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
      }
    | groupTuple( by:0, sort:true ) \
    | concatenate_windows \
    | normalise_vcf \
    | reheader_vcf

    def chromosome_vchks = chromosome_vcfs \
    | summarise_vcf

    concatenate_all(
        (chromosome_vcfs.flatten().collect()),
        (chromosome_vchks.collect()),
        'VARIANTS_UNFILTERED'
    )

}

// Step 1 - Genotyping
process genotype_windows {

    input:
    path crams
    path ploidy_file
    path windows_dir
    each window

    output:
    path("${window}.vcf.gz")

    script:
    """
    if [[ "${window}" == "scaffold"* ]]; then
        bcftools mpileup -d 8000 --ignore-RG -R ${windows_dir}/${window} -a AD,DP,SP -Ou -f ${params.ref} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${window}.vcf.gz
    else
        bcftools mpileup -d 8000 --ignore-RG -r ${window} -a AD,DP,SP -Ou -f ${params.ref} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${window}.vcf.gz
    fi
    """
}

// Step 2 - Concatenate based on key values
process concatenate_windows {

    input:
    tuple val(key), path(window_vcfs_in_key)

    output:
    val(key)
    file("${key}.vcf.gz")
    file("${key}.vcf.gz.csi")

    script:
    """
    key_sorted_windows=\$(echo $window_vcfs_in_key | tr ' ' '\n' | sort -t"-" -k2 -n | tr '\n' ' ')
    bcftools concat --threads ${task.cpus} -n -O z -o ${key}.vcf.gz \${key_sorted_windows}
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 3 - Normalise VCF
process normalise_vcf {

    input:
    val(key)
    file('input.vcf.gz')
    file('input.vcf.gz.csi')
    
    output:
    val(key)
    file("${key}.vcf.gz")
    file("${key}.vcf.gz.csi")

    script:
    """
    # BASIC NORMALISATION
    bcftools norm --threads ${task.cpus} --fasta-ref ${params.ref} -O z -o ${key}_tmp.vcf.gz input.vcf.gz
    bcftools index --threads ${task.cpus} ${key}_tmp.vcf.gz

    # REMOVE SPANNING INDELS AND RE-NORMALISE
    bcftools view --threads ${task.cpus} -V indels -e 'ALT="*" | N_ALT>1' ${key}_tmp.vcf.gz \
    | bcftools norm --threads ${task.cpus} -d exact -O z -o ${key}.vcf.gz
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 4 - Reheader VCF and output to per-chromosome directories
process reheader_vcf {

    publishDir "${params.publish_dir}/chroms/${key}", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    val(key)
    file('input.vcf.gz')
    file('input.vcf.gz.csi')
    
    output:
    tuple file("${key}.vcf.gz"), file("${key}.vcf.gz.csi")

    script:
    """
    bcftools query -l input.vcf.gz | xargs -n 1 basename | awk -F '.' '{print \$1}' > samples.list
    bcftools reheader --threads ${task.cpus} --samples samples.list -o ${key}.vcf.gz input.vcf.gz
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 5 - Get summary stats for each per-chromosome VCF
process summarise_vcf {

    input:
    tuple file(vcf), file(csi)

    output:
    file("${vcf.simpleName}.vchk")

    script:
    """
    bcftools query -l $vcf > samples.list
    bcftools stats --samples-file samples.list $vcf > ${vcf.simpleName}.vchk
    """
}

// Step 6 - Concatenate and output all per-chromosome and per-scaffold VCFs and VCHKs
process concatenate_all {
    
    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    path(collected_vcfs), stageAs: "staged_vcfs/*"
    path(collected_vchks), stageAs: "staged_vchks/*"
    val(collection_name)

    output:
    file("${collection_name}.vcf.gz")
    file("${collection_name}.vcf.gz.csi")
    file("${collection_name}.vchk")

    script:
    """
    # CONCATENATE VCFs (CHECK IF AN ADDITIONAL NORM STEP IS NECESSARY)
    bcftools concat --threads ${task.cpus} -n -O z -o ${collection_name}.vcf.gz staged_vcfs/*.vcf.gz
    bcftools index --threads ${task.cpus} ${collection_name}.vcf.gz

    # CONCATENATE VCHKs
    plot-vcfstats --merge staged_vchks/* > ${collection_name}.vchk

    # CHANGE CONCATENATED VCHK METADATA FOR MULTIQC
    sed -i -e 's/This file was produced by plot-vcfstats/This file was produced by bcftools stats/g' \
           -e 's/\\(.*\\)\\*r[^.]*\\.vcf\\.gz/\\1${collection_name}\\.vcf\\.gz/g' ${collection_name}.vchk 
    """
}
