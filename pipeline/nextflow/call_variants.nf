#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Call VCF
//
// Developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Default input parameters
params.ref = file('/cluster/projects/nn10082k/ref/house_sparrow_genome_assembly-18-11-14_masked.fa')
params.publish_dir = './output'
params.windows_dir = './output/called_variants/genome_windows/'
params.ploidy_file = file('./pipeline/defaults/default.ploidy')

// Workflow
workflow{    

    windows_list = file("${params.windows_dir}/genome_windows.list").readLines()
        
    Channel.fromPath(params.ploidy_file).set{ploidy_file}
    Channel.fromPath(params.windows_dir).set{windows_dir}
    Channel.fromList(windows_list).set{windows}
    Channel.from(file(params.crams)).set{crams}

    genotyping(crams, ploidy_file, windows_dir, windows) \
    | map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
      }
    | groupTuple( by:0,sort:true ) \
    | vcf_concat | vcf_normalise | vcf_reheader | rm_spanning_indels

}

// Step 1 - Genotyping
process genotyping {

    input:
    path crams
    path ploidy_file
    path windows_dir
    each windows

    output:
    path ("${windows}.vcf.gz")

    script:
    """
    if [[ "${windows}" == "scaff"* ]];
    then
        # if window is a scaffold
        bcftools mpileup -d 8000 --ignore-RG -R ${windows_dir}/${windows} -a AD,DP,SP -Ou -f ${params.ref} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${windows}.vcf.gz
    else
        # for normal genome windows
        bcftools mpileup -d 8000 --ignore-RG -r ${windows} -a AD,DP,SP -Ou -f ${params.ref} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${windows}.vcf.gz
    fi
    """
}

// Step 2 - Concatenate based on key values
process vcf_concat {

    input:
    tuple val(key), file(vcfs)

    output:
    tuple \
    val(key), \
    file ("${key}_concat.vcf.gz"), \
    file ("${key}_concat.vcf.gz.csi")

    script:
    """
    # sort the vcfs first 
    sort_vcfs=\$(echo $vcfs | tr ' ' '\n' | sort -t"-" -k2 -n | tr '\n' ' ')
    # then run bctools
    bcftools concat --threads ${task.cpus} -n -O z -o ${key}_concat.vcf.gz \${sort_vcfs}
    bcftools index --threads ${task.cpus} ${key}_concat.vcf.gz
    """
}

// Step 3 - Normalise VCF
process vcf_normalise {

    input:
    tuple \
    val(key),
    file ("${key}_concat.vcf.gz"), \
    file ("${key}_concat.vcf.gz.csi")
    
    output:
    tuple \
    val(key), \
    file ("${key}_norm.vcf.gz"), \
    file ("${key}_norm.vcf.gz.csi")

    script:
    """
    bcftools norm --threads ${task.cpus} --fasta-ref ${params.ref} -O z -o ${key}_norm.vcf.gz ${key}_concat.vcf.gz
    bcftools index --threads ${task.cpus} ${key}_norm.vcf.gz
    """
}

// Step 4 - Reheader VCF
process vcf_reheader {

    input:
    tuple \
    val(key),
    file ("${key}_norm.vcf.gz"), \
    file ("${key}_norm.vcf.gz.csi")
    
    output:
    tuple file ("${key}.vcf.gz"), file ("${key}.vcf.gz.csi")

    script:
    """
    bcftools query -l ${key}_norm.vcf.gz | xargs -n 1 basename | awk -F '_' '{print \$1}' > samples
    bcftools reheader --threads ${task.cpus} -s samples -o ${key}.vcf.gz ${key}_norm.vcf.gz
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 5 - Remove spanning indels and re-normalise
process rm_spanning_indels {
    
    publishDir "${params.publish_dir}/vcf", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    tuple file(vcf), file(index)

    output:
    tuple file ("${vcf.simpleName}.vcf.gz"), file ("${vcf.simpleName}.vcf.gz.csi")

    script:
    """
    bcftools view --threads ${task.cpus} -V indels -e 'ALT="*" | N_ALT>1' $vcf \
    | bcftools norm --threads ${task.cpus} -D -O z -o ${vcf.simpleName}.vcf.gz
    bcftools index --threads ${task.cpus} ${vcf.simpleName}.vcf.gz
    """
}
