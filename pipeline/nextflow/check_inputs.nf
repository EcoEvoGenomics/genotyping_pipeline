#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Check inputs
//
// Pipeline originally developed by Mark Ravinet
// This workflow developed and maintained by Erik Sandertun Røed

workflow {

    assert_contig_names_are_alphanumeric(params.ref_index)

}

process assert_contig_names_are_alphanumeric {

    cpus { 1 }
    memory { 1.GB }
    time { 15.m }

    input:
    path(ref_index)

    script:
    """
    while read -r contig size pos base_per_line byte_per_line; do
        if grep '^[0-9a-zA-Z]*\$' <<< \$contig;
            then : ;
            else
                echo "ERROR: Chromosome name \$contig contains non-alphanumeric characters."
                exit 1
        fi
    done < ${ref_index}
    """
}
