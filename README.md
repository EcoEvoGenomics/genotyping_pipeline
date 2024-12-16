# Variant calling pipeline

Duplicated from https://github.com/markravinet/genotyping_pipeline_v2.git on Wednesday 11 Dec 2024.

### Mark Ravinet 
### 18/11/2024

## Introduction

Welcome to the readme and guide for the Sparrow Ecological & Evolutionary Genommics group variant calling pipeline. After many years of calling variants using a set of multiple slurm scripts, I decided to write a nextflow pipeline that vastly simplifies the process, meaning that you will only have to submit three scripts to go from raw reads to a final, filtered `vcf` ready for analysis.

The scripts are designed to require minimal input from a user but they will produce finished output, statistics and reports for what programs where used in the process. As well as simplifying the process, these scripts are designed to **standardise** the variant calling used within the group, i.e. to prevent different projects using different pipelines and leading to incompatible datasets. However, more importantly, the emphasis here is to improve **reproducibility** - i.e. promoting open science and making it possible for anyone interested to reproduce results in a straightforward way.

There are three steps in the pipeline, each run using a specific nextflow script:

1. `trim_map_realign` - trim reads for quality and adapter sequences, map to a reference genome and perform indel realignment.
2. `call_variants` - call variants across genome windows, concatenate them together to produce per chromosome vcf files.
3. `filter_variants` - filter vcfs for population structure and genome scan analysis.

In addition to these scripts, there are three helper scripts which are optional. These are:

1. `downsample_bams` - a nextflow script that will calculate the depth of a bam, see if it exceeds a threshold you set and if so, downsample it to match that. This is useful if you are combining samples with different read depths. 
2. `create_genome_windows` - a shell script to split your genome of choice into windows for `call_variants` to run on. So you can split the genome into 10 Mb windows and call variants much more efficiently.
3. `concat_vcfs` - concatenate the variant vcfs created by `filter_variants` into a single vcf for the entire genome

## Installation

The simplest way to install everything you need for these scripts to work is to use template `conda` environment. This means you first need to install and setup `conda`. 

**For users working on the Saga infrastructure** you do not need to install `conda`. You can instead follow the guidelines [here](https://github.com/markravinet/markravinet.github.io/wiki/5_using_conda_Saga) to run the pipeline using the shared installed `conda` modules. 

### Installing conda

The full `conda` installation is not necessary. Instead you can use miniconda. Go here and copy the link address for the latest release for Linux 64-bit. Then use `wget` to download it like so:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

From your home directory, simply run this script:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompts and this will install `conda`. Once you have done that, you will need to update it to ensure there are no issues.

```
conda update conda
```

We are then read to set up the `conda` environment you need to run the scripts.

### Using mamba

Conda can be quite slow, so a way to speed it up is to use `mamba`. This is easy to install. You can do so like this:

```
conda install -c conda-forge mamba
```

Once `mamba` is installed you are ready to configure the `conda` environment to use the pipeline.

### Configuring the conda environment

In order to ensure you have all the tools needed to run the pipeline, you can configure the `conda` environment using the `nextflow.yml` provided here. This is also a reference of all the software versions used to run the code. It is very simple to do this, just use the following command:

```
mamba env create --name nf -f nextflow.yml 
```

This will create a `mamba`/`conda` environment called `nf` that contains all the software necessary to run the scripts with minimal manual configuration. 

Finally in order to ensure `abra2` can access the libraries it needs to run, you will need to add the following line to your `bash_profile` or `bashrc`:

```
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib
```

You might also need to increase the amount of memory that `abra2` has access to. In order to do this, you should use `nano` to edit the `~/miniconda3/envs/nf/share/abra2-2.23-1/abra2` so that the `JAVA_TOOL_OPTIONS` variable is increased like so:

```
JAVA_TOOL_OPTIONS="-Xmx8G"
```

### A quick note on Java based errors

Although most of the tools for the pipeline should run without issue, you might get an error that some of the tools cannot run because there is no `java` environment. If this is the case, you need to either install `java` or load the `java` module on your cluster - the pipeline should then work without issue.

## Step 1  - trimming, mapping and aligning

This first script will take your raw reads and run them through [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for a quality assessment. It will then trim them for low-quality bases and remove any adapter sequences. It will then map them to a reference genome of your choice (default is the 2014 House sparrow reference). It also ensures that the samples are renamed to the correct name for all downstream analysis. Finally, the script will produce statistics on the mapping efficiency and depth of coverage of each mapped individual. Note that as of November 2024, the pipeline has been updated so that it only outputs [cram files](https://en.wikipedia.org/wiki/CRAM_(file_format)). 

### The input csv format

In order to run this script, you need to provide it with a csv file with the following columns:

1. Original sample name (only necessary for older samples, otherwise this can just be the sample name)
2. New sample name (i.e. UoN code) - if this is already the sample name in 1, then just repeat it here
3. Forward read location - this should be the **full path** to the forward read
4. Reverse read location - this should be the **full path** to the reverse read
5. Adapter name - the name for the adapter sequences - this should be one of the set of options outlined below.

**A quick note on sample names!** The script will split names based on underscores, so if your samples are called `FH_99` and `FH_100` for example, it will merge them all into a sample called `FH`. So to avoid this, you need to alter them so that they are called `FH99` and `FH100`

Adapters should be written exactly as one of the following options:
- `Illumina_UD-PE`
- `NexteraPE-PE`
- `TruSeq2-PE`
- `TruSeq3-PE-2`
- `TruSeq3-PE`
- `TruSeqUD-PE`

**Note** If your individual has multiple forward and reverse reads, these should all be entered into a separate row for the csv (i.e. for each forward and reverse pair). This will happen when an individual is sequenced across multiple lanes for example. The script is built to account for this and you can give each line the same sample name - **however** the read locations must be different for each. 

An example of the file format is shown below - note that the true file **should not have headers**:
| old name | new name | read1 path | read2 path | adapter |
|------------------|------------------|-------------------------|-------------------------|------------|
| old_sample1_name | new_sample1_name | /path/to/forward_read1/ | /path/to/reverse_read1/ | TruSeq3-PE |
| old_sample2_name | new_sample2_name | /path/to/forward_read2/ | /path/to/reverse_read2/ | TruSeqUD-PE |

### Running the script

Once you have your input csv file ready, you can run the script. To do this, you just need to do the following:

```
nextflow run 1_trim_map_realign.nf --samples samples_test.csv
```

Note that if you run the script in this way, it will run in default mode and use the house sparrow reference genome. In order to alter that, you can add an additional option, `--ref` and specify the location of this alternative reference. Note that in order to do that, you **must** ensure that reference genome has been indexed by bwa prior to running the analysis. There is info on [how to do that here](https://speciationgenomics.github.io/mapping_reference/):

```
nextflow run 1_trim_map_realign.nf --samples samples_test.csv --ref /path/to/alt/ref_genome.fa
```
### Resuming scripts if there are issues

On occassion your script might run into issues. For tips on how to resolve these, see the **Troubleshooting** section below. Once you have fixed this, you can actually restart the nextflow script from a specific point using the `-resume` command. For example:

```
nextflow run 1_trim_map_realign.nf --samples samples_test.csv -resume
```

Remember the script might not pick up exactly where it left off as it runs from preset checkpoints. However, it will generally be much faster doing this than starting from scratch.

### Script outputs

Once complete, the script will create two directories with outputs in:

1. `align` - this contains the mapped, realigned and sorted cramfiles with their indexes for each individual you ran the script on.
2. `stats` - this directory contains statistics from each of the cramfiles for their mapping success and depth of coverage. More information below.

## Step 1a - Downsampling bams

Occassionally, you might need to combine data from sequencing runs with different depths. If this is the case, you will need to downsample individuals sequenced at a higher depth when combining them all together as otherwise this will bias your analysis. For example, in a PCA with 10X indiviudals vs 30X individuals, you will find differences between them that are driven by this variation in sequencing depth. The script `1a_downsample_bams.nf` is designed to handle this and is very easy to run and use.

All you need to do is provide it with a list of bam files, identical to that you would use for the `2_call_variants.nf` (see **Creating a list of bams** below). For example:

```
nextflow run 1a_downsample_bams.nf --bams bams.list --depth 10
```

You should also provide the `--depth` option. The script will then calculate the depth of each sample, see if it is greater than 10X and if so, it will downsample the individuals to this level of coverage using `samtools`. If you don't add a `--depth` option, the script will default to 10 anyway but you can set this to whatever you wish.

Once it has run, the script will only write outputs for the bams which exceeed the specified depth. These downsampled bams and their indexes will be written to the `align` directory with the suffix `_ds.bam` and `_ds.bam.bai` respectively. You can then use these in your variant calling analysis.

## Step 2 - Variant calling

The second script in the pipeline will take a list of bamfiles and performs genotyping on all individuals. To do this, it uses `bcftools` and will call sites at every position in the genome (i.e. it calls invariant sites as well as variants). This is obviously a large job, especially on larger genomes. So to increase efficiency, the script parallelises across genome windows. The default is 10 Mb but you can set these to whatever size you wish. However, tweaking windows has to be done with a separate bash script (see below), not within the nextflow pipeline. After calling genotypes in windows, the script will take care of sorting and concatenating the windows together so that you are left with a vcf file for each chromosome, the mtDNA and also the unanchored scaffolds in your genome.

**NB** as of November 2024, the pipeline has been adapted to use cramfiles instead of bamfiles as these are more efficient in terms of storage space. However, it is still possible to use the `2_call_variants.nf` script with bams if needed. Just point the script to the files and it will work.

### Setting up genome windows

In order to parallelise across genome windows, the script requires a list of said windows in a text file. These are easy to generate using the helper script `0_create_genome_windows.sh`. This simple bash script uses bedtools to split the genome into windows of whatever size you wish (I recommend 10,000,000 bases for most bird genomes) and then it will also create separate files to account for all the scaffolds.

You need to run this script in the directory you are running the nextflow analysis in. You run it like so:

```
bash 0_create_genome_windows.sh ref_index window_size output_name
```

A more specific example, that will work for the house sparrow genome looks like this:

```
bash 0_create_genome_windows.sh /share/Passer/data/reference/house_sparrow_ref.fa.fai 10000000 sparrow_genome_windows.list
```

So option 1 is the path to the reference index, in `fai` format, option 2 is the window size (10 Mb here) and option 3 is the name of the output.

Running this script will produce a set of different files. The first will be text file with a list of all the windows, i.e. `sparrow_genome_windows.list` in the example above. The second will be 10 files called `scaffolds:00`, `scaffolds:01` and so on. The pipeline needs all these files to be present in the base directory you are running nextflow in. 

**NB this script is only really tested on the house sparrow genome** - if it does not work for your genome, you will need to edit it. 

### Creating a list of crams

The other input this script needs is a list of cramfiles. This is very simple - it is just a list of paths of the files that you intend to analyse. If you generated these using `1_trim_map_realign.nf` then the names should already be standardised as `samplename_realigned.cram`. Provided your crams are named this way, then the calling script will also ensure that the sample names are written into the final vcf. Here is an example of what the file should look like:

```
/path/to/align/PDOMNOR8934547_realigned.cram
/path/to/align/PDOMNOR8L19766_realigned.cram
/path/to/align/PDOMNOR8L19786_realigned.cram
/path/to/align/PDOMNOR8L52141_realigned.cram
/path/to/align/PDOMNOR8L52830_realigned.cram
```

Note that these crams do not need top be in the same directory and they do not need to be in a directory called align. This means you can call a vcf from crams in multiple locations easily. 

### Running the script

Once you have your list of crams and genome windows file ready, you can run the script like this:

```
nextflow run 2_call_variants.nf --bams crams.list --windows sparrow_genome_windows.list
```

As above, you can use the `--ref` option to set the location of a specific reference genome. Furthermore, as with all the scripts, you can use the `-resume` option to rerun from a checkpoint if it fails for any reason.

### Optional - specifying ploidy

In some cases, you might need to specify the ploidy of the chromosomes or scaffolds you are calling in your analysis. A typical example of this is when calling SNPs in the mitochondrial genome (i.e. where ploidy is haploid). The pipeline has now been updated to incorporate this functionality but it is **optional**. This means that if you don't specify a ploidy file to `bcftools` then the pipeline will just assume everything is diploid. 

The basic and easiest way to ensure mtDNA is called as haploid is to provide a file that looks like this in the case of the *Passer domesticus* genome:

```
mtDNA 1 16809 F 1
mtDNA 1 16809 M 1

```

This ploidy file is space delimited with the chromosome identity, the start position, the end position, sex, and the ploidy (where 1 = haploid, 2 = diploid). There is more info (here)[https://samtools.github.io/bcftools/bcftools.html#call] on this format. With the example above, the pipeline will assume all other chromosomes are diploid.

You can then run the pipeline exactly as before but this time with the additional `--ploidyFile` option, e.g:

```
nextflow run 2_call_variants.nf --bams bams.list --windows sparrow_genome_windows.list \
--ploidyFile my_ploidy_file.ploidy
```

### Script outputs

The outputs for this script are much simpler than the previous step - it will create a directory called `vcf` and inside will be the gz compressed vcf files for each chromosome, the mtDNA and the genome unanchored scaffolds. These will be raw (i.e. unfiltered) and will contain calls for all variant and non-variant sites in the genome. The directory will also include the indexes for these vcfs.

## Step 3 - filtering

The final script takes control of filtering your vcf files and prepares them for downstream analysis. First it normalises them to remove any issues from concatenating across windows in the calling step. Then it applies custom filters using `vcftools` to create two sets of vcfs; one for population structure analyses (i.e. variants only) and one for genome scan analyses (i.e. variant and invariant sites). There is more information on the outputs below.

### Running the script

The filtering script is the easiest of the three main pipeline scripts to run. It does not require any input as it will automatically look for any vcfs which are gzipped, indexed and which are stored in a directory called `./vcf` in the base directory it is run in. This means it can simply be run with the default filtering options like so:

```
nextflow run 3_filter_variants.nf 
```

However, [as shown here](https://speciationgenomics.github.io/filtering_vcfs/), it is not a good idea to just run filters without checking whether they apply to your dataset. Instead you are able to tweak the filters with a number of options. These are simply provided to the script using option flags and are modified versions of the options for [vcftools](https://vcftools.github.io/examples.html)

- `--miss` - set the missing data at a value between 0 and 1 (where 0 allows 100% missing data and 1 means no missing data); default is 0.8
- `--q_site1` - site quality threshold (as a phred score) for the population structure vcfs - default is 30
- `--q_site2` - site quality threshold (as a phred score) for the genome scan vcfs - default is 30
- `--min_depth` - minimum mean depth of coverage for a variant across all samples - default is 5
- `--max_depth` - maximum mean depth of coverage for a variant across all samples - default is 30
- `--min_geno_depth` - minimum genotype depth per sample. If lower than this value, the genotype will be converted to a missing site - default is 5
- `--max_geno_depth` - maximum genotype depth per sample. If lower than this value, the genotype will be converted to a missing site - default is 30

You can provide all or some of these options to the script using these options. A fully worked example is below:

```
nextflow run 3_filter_variants.nf --miss 0.5 --q_site1 30 --q_site2 40 --min_depth 5 --max_depth 15 --min_geno_depth 5 --max_geno_depth 15
```

You do not need to provide all the options - for example, if you want to just alter the missing data threshold, the following will work.

```
nextflow run 3_filter_variants.nf --miss 0.5 
```

You can also now filter for individuals. This is optional, if you run `3_filter_variants.nf` without any options, it will just include all individuals in the vcf. However if you use `--keep`, the pipeline return a vcf with only those individuals listed in the input file. An example of how to run this is like so:

```
nextflow run 3_filter_variants.nf --keep /path/to/my_inds.txt
```

Here the file provided to the `--keep` option is the same as that used by [`vcftools`](https://vcftools.github.io/man_latest.html) - i.e. this should be a list of individuals like so:

```
Ind1
Ind2
Ind3
```

One thing that is very important if you use the `--keep` option. You **must** provide a full path to the file - i.e. `/path/to/my_inds.txt` - otherwise it will not filter for these individuals.

### Script outputs

The script will create a directory called `vcf_filtered`. Inside this vcf will be per chromosome (and scaffold) vcfs with two different suffixes.

- `chrXX_norm_filtered_ps.vcf.gz` - this is the population structure analysis vcf - it contains variant biallelic SNPs only - ready for PCA, ADMIXTURE and so on.
- `chrXX_norm_filtered_gs.vcf.gz` - this is the genome scan vcf - it contains variant and invariant sites - it is ready for both selection and introgression scans

### Creating a variants whole genome vcf

The `3_filter_variants.nf` script will produce outputs per chromosome. This is done for downstream efficiency but it is sometimes useful to have a whole genome vcf for variants - i.e. the `chrXX_norm_filtered_ps.vcf.gz` vcfs. To produce this, you can use the script `4_concat_vcfs.slurm` which will concatenate all the chromosomes together and normalise the output. This is a slurm script that can be submitted but it needs editing first. You should edit the three variables declared at the start of the script:

- `VCF_FILE` - a text file with the paths of the vcfs to concatenate into a single vcf - this MUST be in order
- `OUTPUT_VCF1` - the name of the output vcf prior to normalisation
- `OUTPUT_VCF2` - the name of the output vcf following normalisation

The `VCF_FILE` should look like this:

```
/path/to/my/vcfs/chr1_norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr2_norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr3_norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr4_norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr5_norm_filtered_ps.vcf.gz
```

The order is important here and should match the order of the chromosomes in the reference genome. 

## Submitting the nextflow pipeline via slurm

The nextflow scripts will take care of the individual scheduling of jobs, but you will need to submit a management job to the cluster using a standard slurm script. This script should be set to run with a relatively low memory (i.e. 4-8Gb), a single CPU and a relatively long time for the entire pipeline to run. For example, you could run it for a week. 

The actual control for each of the different jobs that the nextflow pipeline submits, you should look at the `nextflow.config` file. This is already to set up to use slurm. However, slurm job schedulers will differ between institutions. For the Nottingham Augusta HPC, the `nextflow.config` file will work straight out of the box as it inherits the queue and settings you submit in the master script. For SAGA, you need to edit the `nextflow.config` file to get it to work. You will need to add the `account` value for all of the processes like so, where `XXX` is the account you use on the cluster.

```
   withName: trimming{
   clusterOptions = " --account=XXX --job-name=trim --time=12:00:00 --mem-per-cpu=20G --cpus-per-task=1"
   }
```

**NB if you want to test the pipeline locally** you need to change the name of the `nextflow.config` file so that it is not read by the pipeline automatically. For example, just do this:

```
mv nextflow.config nextflow_config
```

Just remember that when you do go back to submitting the nextflow script, you must change this back so that the pipeline will resume the slurm execution. 

**NB** For group members and users on Saga, you can use the `run_genotype_pipeline_saga.slurm` as a guideline. I have included it on this repo as a template for others too. 

## Generating a quality control report

One aspect of the pipeline is that it is designed to produce outputs that can be easily summarised in a simple quality control report. For this, I recommend using `multiqc`. Once the pipeline is finished (or at any of the intermediate main three steps). 

All you need to do is be in the directory you have run the genotyping pipeline in and then run the following:

```
multiqc .
```

