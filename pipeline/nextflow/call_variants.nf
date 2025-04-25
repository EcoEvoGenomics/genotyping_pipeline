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

    // Input files
    def ref_ploidy = file(params.ref_ploidy_file)
    def ref_genome = file(params.ref_genome)
    def ref_index = file(params.ref_index)

    // Prepare Channels for the genome windows
    define_windows(ref_index, params.window_size, params.ref_scaffold_name)
    def windows_dir = define_windows.out.windows_dir
    def window_list = define_windows.out.windows.map{path -> file(path.toString())}.readLines()
    
    // Prepare a Channel of input CRAM paths
    def crams = Channel.fromPath("${params.cram_dir}/**.cram")
    .map { cram -> "${cram.toString()}\n" }
    .collectFile()

    // Genotype CRAMs in windows and concatenate windows into per-chromsome VCFs
    def chromosome_vcfs = genotype_windows(crams, ref_ploidy, ref_genome, ref_index, windows_dir, window_list) \
    | map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
    } \
    | groupTuple( by:0, sort:true ) \
    | concatenate_windows

    // Normalise and reheader chromosome VCFs
    chromosome_vcfs = normalise_vcf(chromosome_vcfs, ref_genome, ref_index) \
    | reheader_vcf

    // Run bcftools stats on each VCF 
    def chromosome_vchks = chromosome_vcfs \
    | summarise_vcf

    // Concatenate and output VCHKs
    concatenate_vchks(chromosome_vchks.collect(), 'variants_unfiltered')

    // If concatenated VCF wanted, output
    if (params.concatenate_vcf == 'yes') {
        concatenate_vcfs(chromosome_vcfs.flatten().collect(), 'variants_unfiltered')
    }
}

// Step 0 - Write genotyping genome windows
process define_windows {

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

// Step 1 - Genotyping
process genotype_windows {

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
    if [[ "${window}" == "scaffold"* ]]; then
        bcftools mpileup --threads ${task.cpus} -d 8000 --ignore-RG -R ./genome_windows/${window} -a AD,DP,SP -Ou -f ${ref_genome} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${window}.vcf.gz
    else
        bcftools mpileup --threads ${task.cpus} -d 8000 --ignore-RG -r ${window} -a AD,DP,SP -Ou -f ${ref_genome} -b ${crams} \
        | bcftools call --threads ${task.cpus} --ploidy-file ${ploidy_file} -f GQ,GP -mO z -o ${window}.vcf.gz
    fi
    """
}

// Step 2 - Concatenate based on key values
process concatenate_windows {

    input:
    tuple val(key), path(window_vcfs_in_key)

    output:
    tuple val(key), file("${key}.vcf.gz"), file("${key}.vcf.gz.csi")

    script:
    """
    key_sorted_windows=\$(echo ${window_vcfs_in_key} | tr ' ' '\n' | sort -t"-" -k2 -n | tr '\n' ' ')
    bcftools concat --threads ${task.cpus} -n -O z -o ${key}.vcf.gz \${key_sorted_windows}
    bcftools index --threads ${task.cpus} ${key}.vcf.gz
    """
}

// Step 3 - Normalise VCF
process normalise_vcf {

    input:
    tuple val(key), file('input.vcf.gz'), file('input.vcf.gz.csi')
    path(ref_genome)
    path(ref_index)
    
    output:
    tuple val(key), file("${key}.vcf.gz"), file("${key}.vcf.gz.csi")

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

    input:
    tuple val(key), file('input.vcf.gz'), file('input.vcf.gz.csi')
    
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
    bcftools query -l ${vcf} > samples.list
    bcftools stats --samples-file samples.list ${vcf} > ${vcf.simpleName}.vchk
    """
}

// Step 6A - Collect VCHKs and output combined file for MultiQC
process concatenate_vchks {
    
    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    path(collected_vchks), stageAs: "staged_vchks/*"
    val(collection_name)

    output:
    file("${collection_name}.vchk")

    script:
    """
    # CONCATENATE VCHKs
    plot-vcfstats --merge staged_vchks/* > ${collection_name}.vchk

    # CHANGE CONCATENATED VCHK METADATA FOR MULTIQC
    sed -i -e 's/This file was produced by plot-vcfstats/This file was produced by bcftools stats/g' \
           -e 's/\\(.*\\)\\*r[^.]*\\.vcf\\.gz/\\1${collection_name}\\.vcf\\.gz/g' ${collection_name}.vchk 
    """
}

// Step 6B (Optional) - Collect VCFs and output combined file
process concatenate_vcfs {
    
    publishDir "${params.publish_dir}", saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    path(collected_vcfs), stageAs: "staged_vcfs/*"
    val(collection_name)

    output:
    file("${collection_name}.vcf.gz")
    file("${collection_name}.vcf.gz.csi")

    script:
    """
    # CONCATENATE VCFs (CHECK IF AN ADDITIONAL NORM STEP IS NECESSARY)
    bcftools concat --threads ${task.cpus} -n -O z -o ${collection_name}.vcf.gz staged_vcfs/*.vcf.gz
    bcftools index --threads ${task.cpus} ${collection_name}.vcf.gz
    """
}
