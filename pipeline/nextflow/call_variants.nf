#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Call VCF
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

// Workflow
workflow{    

    // Reference files are passed as parameters
    def ref_ploidy = file(params.ref_ploidy_file)
    def ref_genome = file(params.ref_genome)
    def ref_index = file(params.ref_index)
    
    // The CRAMs are parsed to variant calling as a sorted list of paths
    def input_crams = Channel.fromPath("${params.cram_dir}/**.cram").toList()
    def cram_list_file = input_crams \
    | flatMap { cram -> cram } \
    | map { cram -> cram.toString() + "\n" } \
    | collectFile()
    sorted_cram_list_file = sort_cramlist(cram_list_file)

    // The number of CRAMs determines size / number of genotyping windows
    def genotyping_window_base_size = 10000000
    def cram_bucketsize = 75
    def cram_count = input_crams.map { crams -> crams.size() as Integer }
    def cram_bucketcount = cram_count.map { n -> (n + cram_bucketsize - 1) / cram_bucketsize }
    def window_size = cram_bucketcount.map { bc -> (genotyping_window_base_size / bc) as Integer }

    // Genome windows are scaled down from a base size given number of CRAMs
    define_windows(ref_index, window_size, params.ref_scaffold_name)
    def windows_dir = define_windows.out.windows_dir
    def window_list = define_windows.out.windows.map{path -> file(path.toString())}.readLines()

    // Genotyping windows are concatenated into chromosome VCFs
    def chromosome_vcfs = genotype_window(sorted_cram_list_file, ref_ploidy, ref_genome, ref_index, windows_dir, window_list) \
    | map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
    } \
    | groupTuple( by:0, sort:true ) \
    | concatenate_windows

    // Chromosome VCFs are normalised and reheadered
    chromosome_vcfs = normalise_vcf(chromosome_vcfs, ref_genome, ref_index) \
    | reheader_vcf

    // Stats are required per chromosome
    def chromosome_vchks = chromosome_vcfs \
    | summarise_vcf

    // Stats per chromosome are concatenated together
    concatenate_vchks(chromosome_vchks.collect(), "variants_unfiltered")

    // A concatenated VCF is produced if specified in parameters
    if (params.concatenate_vcf == "yes") {
        concatenate_vcfs(chromosome_vcfs.flatten().collect(), ref_index, "", params.ref_scaffold_name, "variants_unfiltered")
    }
}

// Step 0 - Write genotyping genome windows
process define_windows {

    container "quay.io/biocontainers/bedtools:2.30.0--h468198e_3"
    cpus 1
    memory 256.MB
    time 5.m

    input:
    path(ref_index)
    val(window_size)
    val(ref_scaffold_name)

    output:
    path('genome_windows/*'), emit: windows_dir
    path('genome_windows/genome_windows.list'), emit: windows

    script:
    """
    # DIRECTORY OF GENOME WINDOW FILES
    mkdir genome_windows

    # MAKE FILE OF CHROMOSOME AND SCAFFOLD SIZES
    cut -f 1-2 ${ref_index} > genome_windows/genome_size.txt

    # MAKE CHROMOSOME WINDOWS
    cd genome_windows
    bedtools makewindows -g genome_size.txt -w ${window_size} \
    | grep -v ${ref_scaffold_name} \
    | awk '{print \$1":"\$2"-"\$3}' \
    > genome_windows.list

    # MAKE SCAFFOLD WINDOWS
    bedtools makewindows -g genome_size.txt -w ${window_size} \
    | grep ${ref_scaffold_name} \
    | awk '{print \$1":"\$2"-"\$3}' \
    > scaffolds.list

    # FIRST POSITION MUST BE 1, NOT 0
    sed -i 's/:0-/:1-/g' genome_windows.list
    sed -i 's/:0-/:1-/g' scaffolds.list

    # FOR SCAFFOLDS, MAKE TAB-DELIMITED TEMPORARY FILE
    cat scaffolds.list \
    | tr ":" "-" \
    | awk -F "-" '{print \$1"\\t",\$2"\\t",\$3}' \
    > scaffolds.tmp
    mv scaffolds.tmp scaffolds.list

    # CREATE SCAFFOLD GROUP FILES WITH NO MORE THAN window_size BASES TO GENOTYPE
    file_counter=1
    base_counter=0
    while read -r scaffold start end; do

        scaffold_window_bases=\$((end - start))
        
        # If adding window to current file would exceed window size, start new file
        if ((base_counter + scaffold_window_bases > ${window_size})); then
            ((file_counter++))
            base_counter=0
        fi
        
        output_file="scaffolds:\$(printf "%02d" "\$file_counter")"
        echo -e "\$scaffold\\t\$start\\t\$end" >> "\$output_file"

        base_counter=\$((base_counter + scaffold_window_bases))

    done < scaffolds.list

    # ADD SCAFFOLDS SUBGROUPS TO WINDOW LIST ONLY IF ANY EXIST
    n_scaffolds=\$(ls scaffolds:* | wc -l)
    if [ \${n_scaffolds} -gt 0 ]
    then
        for i in scaffolds:*; do echo \$i; done >> genome_windows.list
    fi
    """
}

// Step 0 - Sort CRAM list
process sort_cramlist {

    cpus 1
    memory 256.MB
    time 5.m

    input:
    path(crams)

    output:
    path("crams_sorted.list")

    script:
    """
    sort -V ${crams} > crams_sorted.list
    """
}

// Step 1 - Genotyping
process genotype_window {

    label "require_pipefail"

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus { 1 }
    memory { 4.GB * task.attempt }
    time { 8.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    path(crams)
    path(ploidy_file)
    path(ref_genome)
    path(ref_index)
    path(windows_dir), stageAs:'./genome_windows/'
    each(window)

    output:
    path("${window}.vcf.gz")

    script:
    """
    # CALL NONSCAFFOLD WINDOWS AND SCAFFOLD WINDOWS SEPARATELY
    call_nonscaffold() {
        bcftools mpileup --threads ${task.cpus} -d 8000 --ignore-RG -r ${window} -a AD,DP,SP -Ou -f ${ref_genome} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${window}.vcf.gz
    }
    call_scaffold() {
        bcftools mpileup --threads ${task.cpus} -d 8000 --ignore-RG -R ./genome_windows/${window} -a AD,DP,SP -Ou -f ${ref_genome} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${window}.vcf.gz
    }
    if [[ "${window}" == "scaffold"* ]]; then
        call_scaffold
    else
        call_nonscaffold
    fi
    """
}

// Step 2 - Concatenate based on key values
process concatenate_windows {

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 1
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(key), path(window_vcfs_in_key)

    output:
    tuple val(key), path("${key}.vcf.gz"), path("${key}.vcf.gz.csi")

    script:
    """
    key_sorted_windows=\$(echo ${window_vcfs_in_key} | tr ' ' '\n' | sort -t"-" -k2 -n | tr '\n' ' ')
    bcftools concat --threads ${task.cpus} -n -O z -o ${key}.vcf.gz \${key_sorted_windows}
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 3 - Normalise VCF
process normalise_vcf {

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 16
    memory { 8.GB * task.attempt }
    time { 8.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(key), path('input.vcf.gz'), path('input.vcf.gz.csi')
    path(ref_genome)
    path(ref_index)
    
    output:
    tuple val(key), path("${key}.vcf.gz"), path("${key}.vcf.gz.csi")

    script:
    """
    # BASIC NORMALISATION
    bcftools norm --threads ${task.cpus} --fasta-ref ${ref_genome} -O z -o ${key}_tmp.vcf.gz input.vcf.gz
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

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 1
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple val(key), path('input.vcf.gz'), path('input.vcf.gz.csi')
    
    output:
    tuple path("${key}.vcf.gz"), path("${key}.vcf.gz.csi")

    script:
    """
    bcftools query -l input.vcf.gz | xargs -n 1 basename | awk -F '.' '{print \$1}' > samples.list
    bcftools reheader --threads ${task.cpus} --samples samples.list -o ${key}.vcf.gz input.vcf.gz
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 5A - Get summary stats for each per-chromosome VCF
process summarise_vcf {

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 1
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    tuple path(vcf), path(csi)

    output:
    path("${vcf.simpleName}.vchk")

    script:
    """
    bcftools query -l ${vcf} > samples.list
    bcftools stats --samples-file samples.list ${vcf} > ${vcf.simpleName}.vchk
    """
}

// Step 5B - Collect VCHKs and output combined file for MultiQC
process concatenate_vchks {

    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 1
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3

    input:
    path(collected_vchks), stageAs: "staged_vchks/*"
    val(collection_name)

    output:
    path("${collection_name}.vchk")

    script:
    """
    # CONCATENATE VCHKs
    plot-vcfstats --merge staged_vchks/* > ${collection_name}.vchk

    # CHANGE CONCATENATED VCHK METADATA FOR MULTIQC
    sed -i -e 's/This file was produced by plot-vcfstats/This file was produced by bcftools stats/g' \
           -e 's/\\(.*\\)\\*r[^.]*\\.vcf\\.gz/\\1${collection_name}\\.vcf\\.gz/g' ${collection_name}.vchk 
    """
}

// (Optional) - Collect VCFs and output combined file with duplicates removed
process concatenate_vcfs {
    
    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'

    container "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpus 2
    memory { 1.GB * task.attempt }
    time { 4.h * task.attempt }

    errorStrategy "retry"
    maxRetries 3
    
    input:
    path(collected_vcfs), stageAs: "staged_vcfs/*"
    path(ref_index)
    val(vcf_suffix)
    val(ref_scaffold_name)
    val(collection_name)

    output:
    path("${collection_name}.vcf.gz")
    path("${collection_name}.vcf.gz.csi")

    script:
    """
    cat ${ref_index} | grep -v ${ref_scaffold_name} | awk '{print "./staged_vcfs/" \$1 "${vcf_suffix}.vcf.gz"}' > reference_sorted_vcfs.list
    echo "./staged_vcfs/scaffolds${vcf_suffix}.vcf.gz" >> reference_sorted_vcfs.list
    bcftools concat --threads ${task.cpus} --file-list reference_sorted_vcfs.list --output-type u \
    | bcftools norm --threads ${task.cpus} -d exact --output-type z --output ${collection_name}.vcf.gz
    bcftools index --threads ${task.cpus} ${collection_name}.vcf.gz
    """
}
