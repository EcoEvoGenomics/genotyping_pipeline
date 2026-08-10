#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Check inputs
//
// Pipeline originally developed by Mark Ravinet
// This workflow developed and maintained by Erik Sandertun Røed

workflow {

    def ref_index = file(params.ref_genome.toString() + ".fai")
    check_ref_contig_names(ref_index)

}

process check_ref_contig_names {

    cpus { 1 }
    memory { 1.GB }
    time { 15.m }

    input:
    path(ref_index)

    script:
    """
    while read -r contig size pos base_per_line byte_per_line; do
        if grep -qE '^[[:alnum:]_]+\$' <<< \$contig;
            then : ;
            else
                echo "ERROR: Chromosome name \$contig contains non-alphanumeric, non-underscore characters."
                exit 1
        fi
    done < ${ref_index}
    """
}

