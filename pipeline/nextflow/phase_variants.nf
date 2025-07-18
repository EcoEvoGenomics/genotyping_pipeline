#!/usr/bin/env nextflow

// CEES Ecological and evolutionary genomics group - genotyping pipeline
// https://github.com/EcoEvoGenomics/genotyping_pipeline
//
// Workflow: Phase VCF file
//
// Originally developed by Mark Ravinet
// Co-developed and maintained by Erik Sandertun Røed

include { define_windows } from "./call_variants.nf"

