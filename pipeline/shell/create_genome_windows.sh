#!/bin/sh

# Make genome windows

# Set variables
REF_IDX=$1
WINDOW_SIZE=$2
SCAFFOLD_NAME=$3
OUTPUT_DIR=$4

cd $OUTPUT_DIR

# make genome size file
cut -f 1-2 ${REF_IDX} > genome_size.txt

# make the windows and cat together
bedtools makewindows -g genome_size.txt -w ${WINDOW_SIZE} | \
grep -v ${SCAFFOLD_NAME} | awk '{print $1":"$2"-"$3}' \
> genome_windows.list

# GATK will fail if the coords include 0, so edit to start from 1
sed -i_bak 's/:0-/:1-/g' genome_windows.list

# creates 115 genome windows

# for scaffolds
bedtools makewindows -g genome_size.txt -w ${WINDOW_SIZE} | \
grep ${SCAFFOLD_NAME} | awk '{print $1":"$2"-"$3}' \
> scaffolds.list

# GATK will fail if the coords include 0, so edit to start from 1
sed -i_bak 's/:0-/:1-/g' scaffolds.list

# next need to use awk in order to make this a tab delim file
cat scaffolds.list | tr ":" "-" | awk -F "-" '{print $1"\t",$2"\t",$3}' > scaffolds.list2

# split scaffolds into multiple files
total_lines=$(wc -l <scaffolds.list2)
num_files=10
((lines_per_file = (total_lines + num_files - 1) / num_files))

# Split the actual file, maintaining lines.
split -d --lines=${lines_per_file} scaffolds.list2 scaffolds_

# add scafs to windows list only if any exist
n_scaffolds=$(ls scaffolds_* | wc -l)

if [ ${n_scaffolds} -gt 0 ]
then
    for i in scaffolds_*; do echo $i; done >> genome_windows.list
fi
