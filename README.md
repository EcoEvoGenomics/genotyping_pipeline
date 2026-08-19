# XENO: Expedient Genotyping for Nonmodel Organisms

## What is XENO?
XENO is a beginner-friendly genotyping pipeline built for evolutionary and ecological genomics. With XENO, you can (1) trim reads and align them to a reference genome, (2) call variants, (3) filter those variants, and finally (4) phase the variants. These four stages of XENO may be performed separately in stepwise order or end-to-end. At each stage, XENO outputs alignments, variant call files, and / or quality control statistics. These outputs may be passed directly to workflows from XENO's companion repository [RIPLEY](https://github.com/EcoEvoGenomics/RIPLEY), which automate many of the most common analyses in evolutionary and ecological genomics.

The Wiki describes XENO in greater detail. There you will also find guides for some specific use-cases or users.

## Installation
### Prerequisites
*In short: Nextflow 25.04.6, Slurm, Apptainer, Internet access.*

XENO is intended to run on a high-performance computing (HPC) system and is built with [Nextflow](https://www.nextflow.io/) version 25.04.6. This *specific* version of Nextflow must be installed to use XENO. Out of the box, XENO also depends on the HPC job manager Slurm and the software container system Apptainer. Users familiar with Nextflow may alternatively reconfigure the nextflow.config file to use another job manager and another Docker-compatible container system, but regrettably we are unable to provide much support for such use cases. Once you have installed Nextflow on an HPC with Slurm and Apptainer, XENO automatically fetches the remaining software requirements. Note that XENO therefore must be run from an environment with internet access.

### Installing Nextflow with Conda
You may install Nextflow 25.04.6 any way you prefer to run XENO. If you wish, you can replicate the Conda environment we used to develop XENO with the included [YAML file](https://github.com/EcoEvoGenomics/XENO/blob/main/examples/nextflow-25.04.6.conda.yaml). You may copy that file directly from GitHub or run the following command after [installing XENO](#installing-xeno).
```sh
conda create --name nf --file examples/nextflow-25.04.6.conda.yaml 
```

### Installing XENO
The intended way to "install" XENO is simply to clone this GitHub repository to a suitable location on your HPC environment:
```sh
git clone https://github.com/EcoEvoGenomics/XENO
```
Be mindful that XENO can and often will produce terabytes of working files. It is therefore ill-advised to install XENO to a location with limited storage space and especially to locations where you share limited storage space with others (unless you are certain what you are doing).

### Installing a development version of XENO
While XENO is not in active development for new features, we do strive to maintain it. For instance we try to respond to bug reports and improve robustness to different inputs as more users bring their unique data to XENO. As such, there will on occasion be branches that are more recent than `main`. If you wish to use an in-development version you can obtain a specific branch, e.g. `dev`:
```sh
git clone -b dev https://github.com/EcoEvoGenomics/XENO
```

## Quickstart
### Summary

1. Ensure XENO has direct read-access to your (`fastq.gz`) sequencing files. If they are stored externally, you must usually download them.
2. Ensure XENO has direct read-access to your reference genome files and that the reference files are correctly formatted (see [here](#required-reference-files)).
3. Prepare an input file with sample information in comma-separated (`.csv`) format (see [here](#how-to-correctly-format-your-input-csv-file)).
4. Adjust options as necessary and finally launch XENO (see [here](#how-to-launch-xeno)).

### Required reference files
XENO requires three mandatory and one optional piece of reference information. They are:

- A reference genome in uncompressed fasta format (`.fa`, `.fasta`). This must have been separately indexed with `bwa index`, and the index output files (`.amb`, `.ann`, `.bwt`, `.fai`, `.pac`, and `.sa`) must be found at the same path. Contigs and scaffolds in the reference genome **must** be named with alphanumerical characters and underscores only. For reference genomes that do not conform to this requirement, simply rename the contigs in the reference index and in supporting metadata such as GFF files. We regret the inconvenience and may update XENO in the future to accomodate more contig names.
- A prefix (as a text string) to distinguish scaffolds from contigs. Many reference genomes distinguish loosely assembled scaffolds from full contigs (e.g. `chr1`, `chr2`, ...) by prefixing their names with different strings. For instance, NCBI reference genomes may prefix contigs with "NC_" and scaffolds with "NW_". If applicable, you should provide the prefix characteristic of scaffolds (e.g. `NW_`) to enable more efficient distribution of genotyping windows. If not applicable, you may provide an arbitrary string that matches no contig.
- A ploidy file (such as in the `examples/` directory). This file is required to call sex chromosomes, haploid chromosomes, mtDNA, or other non-diploid chromosomes correctly. By default every contig is considered diploid.
- **Optional:** For variant phasing, you *may* provide XENO the absolute path to a directory containing recombination rate maps. The recombination maps must be compatible with [SHAPEIT5](https://odelaneau.github.io/shapeit/). The directory should contain one recombination map file for each contig in the reference genome and they must be named on the format "contig.map". If for instance your reference genome has the contigs "chr1" and "chr2", the recombination map directory must contain two files: "chr1.map" and "chr2.map". If you do not have recombination maps, you may set this variable to an arbitrary path: XENO will statistically phase the variants instead of using maps.

### How to correctly format your input CSV file
For each sample you wish to genotype, XENO requires five inputs. You must provide the inputs in a comma-separated file (`.csv`) with one row per pair of forward and reverse read files and five columns:

1. Sample ID. A single text string of *only* alphanumeric characters, e.g. "SAMPLE1" but not "SAMPLE_1" or "SAMPLE.1".
2. Sample sex. The sex codes (e.g. F, M) are arbitrary but must correspond to the reference ploidy file (see [above](#required-reference-files)).
3. Lane code to distinguish reads from the same sample from different sequencing lanes. Use the format "LXXX" where "XXX" is a number with leading zeroes (e.g. L001). If you have only one set of files per sample, use only "L001". If you have more sets of file per sample, separate them on different rows with increasing lane codes. In *principle* the format (LXXX) is arbitrary, but the QC report will only be organised correctly when this format is adhered to.
4. Absolute path to the forward read file (R1).
5. Absolute path to the reverse read file (R2).

These inputs must be in the above order from left to right. For example:

| Sample ID | Sex | Lane | Path to R1 FASTQ.GZ file | Path to R2 FASTQ.GZ file |
|------------------|---|------|-------------------------|-------------------------|
| PDOM2024IND0001M | M | L001 | /path/to/1M_L001_R1.fastq.gz | /path/to/1M_L001_R2.fastq.gz |
| PDOM2024IND0001M | M | L002 | /path/to/1M_L002_R1.fastq.gz | /path/to/1M_L002_R2.fastq.gz |
| PDOM2024IND0002F | F | L001 | /path/to/2F_R1.fastq.gz | /path/to/2F_R2.fastq.gz |

**Do not include headers** and ensure that the file corresponds to regular CSV conventions: separate values only by commas, add no trailing or leading spaces, and include a newline at the end. You should prepare the CSV in a raw text editor to avoid unexpected formatting discrepancies from spreadsheet software such as Microsoft Excel or Numbers.

### How to launch XENO
#### Configuring Nextflow
Before launching XENO, you must configure the file `nextflow.config` so Nextflow interfaces correctly with the HPC system you use. For instance, you should configure Nextflow to use the correct job queue names for jobs which require access to GPUs or large amounts of memory:
```yaml
process {

  withLabel: "require_gpu" {
    clusterOptions = "--account=your_account --gpus=1"
    queue = "a100"
    cpus = 8
  }

  withLabel: "high_mem_per_cpu" {
    queue = "bigmem"
  }
  
}
```
Consult the documentation for the HPC you use to select an appropriate queue for each of these labels. This is *especially* salient for the `high_mem_per_cpu` label. The cost of jobs with a high memory allocation could increase drastically on a queue with limited memory - at no benefit to you (or other users of the same HPC). While to some extent we have optimised XENO's resource consumption, XENO does not (and cannot) estimate the monetary cost of a job: your use of paid (and indeed, shared) resources is your own responsibility. After configuring Nextflow appropriately, you may optionally modify settings such as the maximum number of concurrent tasks. See the [Nextflow documentation](https://docs.seqera.io/nextflow/config) for options.

#### Setting XENO options
XENO has user-configurable options. For instance you may change the alignment method or variant filters you apply. These options are read from `options.yaml`, which you must configure before launching XENO. The available options are:

| Option | Description | Example | Default |
|--------|-------------|---------|---------|
| `samples` | Path to sample CSV file. | `/user/path/samples.csv` | |
| `trim_align` | Run trimming and alignment stage? | `true` | `true` |
| `call_variants` | Run variant calling stage? | `true` | `true` |
| `filter_variants` | Run variant filtering stage? | `true` | `true` |
| `phase_variants` | Run variant phasing stage? | `true` | `true` |
| `deduplicate` | Deduplicate reads before trimming? | `false` | `false` |
| `downsample` | Downsample R1 and R2 files to `read_target` before trimming? | `false` | `false` |
| `read_target` | Number of reads to downsample to in each of R1 and R2 files.  | `250000` | `1000000` |
| `aligner` | Align with `gpu` ([fq2bam](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam)), `mem` (bwa mem), or `aln` (bwa aln)? While `gpu` is most efficient, you need [compatible GPUs](https://docs.nvidia.com/clara/parabricks/get-started/installation-requirements#hardware-requirements) to use it.  | `mem` | `gpu` |
| `exclude_flags` | Exclude reads with this/these flag(s) from alignments. See options [here](https://www.htslib.org/doc/samtools-flags.html). | `DUP,UNMAP` | `0x400` |
| `concatenate_raw_vcf` | If `false`, only output variants in per-chromosome files. If `true`, also create whole-genome VCF of raw variants. | `true` | `false` |
| `filtering_label` | A label for the filters in `filtering_flags`. Change between runs to re-filter with different settings. | `biallelic_variants` | `default_filters` |
| `filtering_flags` | Path to a file with [VCFtools filtering flags](https://vcftools.github.io/man_latest.html#SITE%20FILTERING%20OPTIONS). Regardless, XENO only retains SNPs - not indels. | `/user/path/filters_biallelic_variants.txt` | `./example/default_filters.txt` |
| `phasing_window_size` | Number of phasing windows. Greater numbers yield greater parallelisation.  | `20000000` | `10000000` |
| `ref_genome` | Path to [reference genome](#required-reference-files). | `/user/path/ref/reference_genome.fa` | |
| `ref_recombination_map_dir` | Path to directory of [recombination maps](#required-reference-files). | `/user/path/ref/recombination_maps/` | |
| `ref_scaffold_name` | [Prefix](#required-reference-files) characteristic of scaffolds in reference genome. | `NW_` | |
| `ref_ploidy_file` | Path to [ploidy file](#required-reference-files). | `/user/path/ref/reference.ploidy` | `./examples/default.ploidy` |

#### Launching XENO through the job scheduler
Experienced UNIX / HPC users will require no assistance launching XENO: simply load the required dependencies and call `bash ./XENO` in the way you see fit. This section and the next will cover two such ways for less experienced bioinformaticians.

A simple way to run XENO is to write a Slurm job script to load dependencies and launch XENO. For example, see `./examples/quickstart.sh`:
```sh
#!/bin/bash

# ADMIN
#SBATCH --job-name=XENO
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k

# RESOURCE ALLOCATION
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=99:00:00

# This quickstart script works on the NRIS Saga HPC for users with access
# to the nn10082k project number. Must be run from top-level in the repo.

module --quiet purge
module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/Nextflow25.04.6
bash ./XENO
```
Note that the `quickstart.sh` script was written for the NRIS Saga HPC. It is a working example, but not universally applicable. You should consult your HPC documentation to modify this script to the specifications of your HPC environment. When you have prepared your Slurm script, simply delegate it to a compute node with:
```sh
sbatch your_launch_script.sh
```
You may then monitor XENO's progress by perusing the slurm logs, and by inspecting the Slurm queue with `squeue --me`. Note that launching XENO through `sbatch` means launching Nextflow on a compute node: this may lead to unexpected problems. If you experience errors while launching XENO through the job scheduler, read on below.

#### Launching XENO on a screen terminal
The most robust way to launch XENO is directly from the terminal. However, XENO should not usually be launched from a regular terminal window, because processes in a regular terminal window are interrupted when you exit. Instead you can open a screen terminal, which persists indefinitely in the background. If you are unfamiliar, the following is a brief introduction. You open a screen terminal with:
```sh
screen
```
You may exit the screen terminal by simultaneously holding the keys Ctrl + A, then the key D on your keyboard. Then, to resume the screen session, use:
```sh
screen -r
```
If you have multiple screen sessions running, you can use the following command to list the active screen terminals:
```sh
screen -list
```
Then you may use the following to resume a specific screen session:
```sh
screen -r session-id-here
```
You must load the necessary prerequesite software on the screen terminal before launching XENO. Like [above](#launching-xeno-through-the-job-scheduler), the exact steps required will depend on your HPC, but this should look something like the following. Please note that these commands will not (necessarily) work - they are for illustration. Consult your HPC documentation to get information about available modules and how to activate Conda.
```sh
module load conda
conda activate Nextflow-25.04.6
bash ./XENO
```
Alternatively, you may [prepare a script with the necessary commands](#launching-xeno-through-the-job-scheduler) and run it directly from the screen terminal (you may omit Slurm headers as they are ignored here):
```sh
bash your_launch_script.sh
```
Either way, on the screen terminal you can monitor Nextflow in real-time as it launches compute jobs. The individual compute jobs are visible through `squeue --me` as usual. The queue is accessible from any terminal. No heavy computation is performed on the screen terminal, which primarily handles software downloads and monitors the compute jobs.

#### Resuming XENO
It is not entirely unusual for HPC jobs to stop prematurely - perhaps because your compute allocation was depleted, because you accidentally cancelled a job, or for altogether enigmatic and irreproducible reasons. Thankfully, Nextflow caches progress (in a subdirectory called `/work`) so that XENO can resume in the event of an untimely stop. Simply ensure all XENO jobs were cancelled, then relaunch. Of course, XENO most often stops in the event of (user) error. XENO is most prone to errors caused by faulty read files for individual samples. If this befalls you, you may simply remove the offending sample from the input CSV and relaunch XENO. We are unable to troubleshoot read files for you, but if XENO persistently exits with an error for another reason, feel free to post a GitHub Issue. Please describe the error in sufficient detail to reproduce it.

## How to cite XENO
Thank you for using XENO. If you wish to cite XENO, you should first cite the third-party software XENO uses. Below is an exhaustive list: please take care to cite appropriately to your own use-case.

### Citing third-party software
- Read deduplication and downsampling: [SeqKit v. 2.10.0](https://doi.org/10.1002/imt2.191)
- Read trimming: [fastp v. 0.24.0](https://doi.org/10.1002/imt2.70078); [FastQC v. 0.12.1](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Alignment ("gpu" option): [Parabricks v. 4.5.0-1](https://doi.org/10.1101/2025.07.23.666378); [SAMtools v. 1.17](https://doi.org/10.1093/gigascience/giab008)
- Alignment ("[mem](https://doi.org/10.48550/arXiv.1303.3997)" and "[aln](https://pubmed.ncbi.nlm.nih.gov/19451168/)" option): [bwa v. 0.7.17](https://doi.org/10.48550/arXiv.1303.3997); [SAMtools v. 1.17](https://doi.org/10.1093/gigascience/giab008); [GATK4 v. 4.6.2.0](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/)
- Variant calling: [BEDtools 2.30.0](https://doi.org/10.1093/bioinformatics/btq033); [BCFtools 1.17](https://doi.org/10.1093/gigascience/giab008)
- Variant filtering: [VCFtools 0.1.16](https://doi.org/10.1093/bioinformatics/btr330); [BCFtools 1.17](https://doi.org/10.1093/gigascience/giab008)
- Variant phasing: [BEDtools 2.30.0](https://doi.org/10.1093/bioinformatics/btq033); [SHAPEIT5 v. 5.1.1](https://doi.org/10.1038/s41588-023-01415-w); [BCFtools 1.17](https://doi.org/10.1093/gigascience/giab008)
- Quality control report: [MultiQC v. 1.28](http://dx.doi.org/10.1093/bioinformatics/btw354)

### Citing XENO
We plan to provide a citeable persistent identifier for XENO later. For the time being, you are welcome to cite this GitHub repository.

______
XENO v. 4.0.0 | 2026 | Erik Sandertun Røed & Mark Ravinet | https://github.com/EcoEvoGenomics/XENO 
