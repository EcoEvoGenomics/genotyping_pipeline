#!/usr/bin/env nextflow

// nextflow pipeline for variant calling
// developed by Mark Ravinet - 08/03/2023
// v 0.3 - 20/06/2023

// script paramaters
params.ref = file('/share/Passer/data/reference/house_sparrow_ref.fa')

// read in a file of bams
//params.bams = file('test_bams.list')
// read in a file of genome windows
// params.windows = file('sparrow_genome_windows.list')

windows_list = file(params.windows)
    .readLines()
    //.each { println it }

// create bams channel
Channel
    .from( file(params.bams))
    .set{bams}
// create windows channel
Channel
    .fromList( windows_list )
    .set{windows}

// ploidy file for genotyping
if (!params.ploidyFile){
    println "No ploidy file provided - by default, all chromosomes will be called as diploid."
    ploidyFile = file('default.ploidy')
} else {
    ploidyFile = file(params.ploidyFile)
}

// View results (to test)
// bams.view()
// windows.view()
// chrs.view()
// result.bams.view()
// result.samples.view()

// genotyping
process genotyping {

    input:
    path (bams)
    path ploidyFile
    each windows

    output:
    path ("${windows}.vcf.gz")

    """
    if [[ "${windows}" == "scaff"* ]];
    then
        # if window is a scaffold
        bcftools mpileup -d 8000 --ignore-RG -R ${baseDir}/${windows} -a AD,DP,SP -Ou -f ${params.ref} -b ${bams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidyFile} -f GQ,GP -mO z -o ${windows}.vcf.gz
    else
        # for normal genome windows
        bcftools mpileup -d 8000 --ignore-RG -r ${windows} -a AD,DP,SP -Ou -f ${params.ref} -b ${bams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidyFile} -f GQ,GP -mO z -o ${windows}.vcf.gz
    fi
    """
}

// concat - based on key value
process vcf_concat {

    input:
    tuple val(key), file(vcfs)

    output:
    tuple \
    val(key), \
    file ("${key}_concat.vcf.gz"), \
    file ("${key}_concat.vcf.gz.csi")

    """
    # sort the vcfs first 
    sort_vcfs=\$(echo $vcfs | tr ' ' '\n' | sort -t"-" -k2 -n | tr '\n' ' ')
    # then run bctools
    bcftools concat --threads ${task.cpus} -n -O z -o ${key}_concat.vcf.gz \${sort_vcfs}
    bcftools index --threads ${task.cpus} ${key}_concat.vcf.gz
    """
}

// normalise vcf
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

    """
    bcftools norm --threads ${task.cpus} --fasta-ref ${params.ref} -O z -o ${key}_norm.vcf.gz ${key}_concat.vcf.gz
    bcftools index --threads ${task.cpus} ${key}_norm.vcf.gz
    """
    
}

// reheader vcf
process vcf_reheader {

    publishDir 'vcf', saveAs: { filename -> "$filename" }

    input:
    tuple \
    val(key),
    file ("${key}_norm.vcf.gz"), \
    file ("${key}_norm.vcf.gz.csi")
    
    output:
    tuple \
    val(key), \
    file ("${key}.vcf.gz"), \
    file ("${key}.vcf.gz.csi")

    """
    bcftools query -l ${key}_norm.vcf.gz | xargs -n 1 basename | awk -F '_' '{print \$1}' > samples
    bcftools reheader --threads ${task.cpus} -s samples -o ${key}.vcf.gz ${key}_norm.vcf.gz
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
    
}

// workflow starts here!

workflow{    
    // set the reference genome from the command line:
    params.ref = "--ref"
    genotyping(bams, ploidyFile, windows) \
    | map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
      }
    | groupTuple( by:0,sort:true ) \
    | vcf_concat | vcf_normalise | vcf_reheader
}
