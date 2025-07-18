#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Phase VCF file
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

include { define_windows } from "./call_variants.nf"

params.vcf = './vcf/sparrows_variants_norm.vcf.gz'
params.map_path = '/cluster/projects/nn10082k/recombination_maps'

// test channel
Channel
    .fromPath("${params.vcf}*") 
    .collect()
    .set{unphased_vcf}

unphased_vcf.view()

// set up windows
windows_list = file(params.windows)
    .readLines()
    //.each { println it }

// make windows channel
Channel
    .fromList( windows_list )
    .set{windows}

// windows.view()
//windows.region.view()
//windows.window.view()
//unphased_vcf.view()

// phase common variants
process phase_common {
    
    errorStrategy 'ignore'
    publishDir 'vcf_phase_window', saveAs: { filename -> "$filename" }

    input:
    path (unphased_vcf)
    each window

    output:
    //path "${window}_phased.bcf"
    tuple \
      path ("${window}_phased.bcf"), \
      path ("${window}_phased.bcf.csi") 

    """
    # create map parameter
    WINDOW=${window}
    echo "This is \${WINDOW} on \${WINDOW%:*}"
    MAP="${params.map_path}/\${WINDOW%:*}.map"

    if [[ -f \${MAP} ]]; then

        echo "\$MAP is present"
        SHAPEIT5_phase_common --input ${unphased_vcf} --region ${window} --map \${MAP} --output ${window}_phased.bcf --thread 12

    else

        echo "\$MAP is not present"
        SHAPEIT5_phase_common --input ${unphased_vcf} --region ${window} --output ${window}_phased.bcf --thread 12

    fi

    """
}

// workflow starts here!
workflow{
    phase_common(unphased_vcf, windows)
}

Channel
    .fromPath( './vcf_phase_window/*.bcf*' )
    .map { file ->
        def key = file.baseName.toString().tokenize(':').get(0)
        return tuple(key, file)
      }
    .groupTuple( by:0,sort:true )
    .set{phased_vcf_windows}

phased_vcf_windows.view()

// ligate phased vcfs - nb this must be done WITHIN chromosomes - it cannot be done across the entire genome
process ligate_chr {

    // publish simlinks into a final vcf directory
    publishDir 'vcf_phased', saveAs: { filename -> "$filename" }, mode: 'copy'

    input:
    tuple val(key), path(vcfs)

    output:
    tuple \
      path ("${key}_phased.bcf"), \
      path ("${key}_phased.bcf.csi") 

    """
    # sort vcfs (and remove indexes)
    sort_vcfs=\$(echo ${vcfs} | tr ' ' '\n' |  grep -v ".csi" | sort -t"-" -k2 -n | tr '\n' ' ')

    # create ligate file
    for i in \$sort_vcfs; do echo \$i >> chunks.txt ; done

    # ligate and produce bcf
    SHAPEIT5_ligate --input chunks.txt --output ${key}_phased.bcf --index --thread 12
    """
}

// workflow starts here!
workflow{
    ligate_chr(phased_vcf_windows) | view
}
