#!/bin/sh

# make genome windows

# set variables
REF_IDX=$1
WINDOW=$2
SCAFFOLD_NAME=$3
OUTPUT=$4
# REF_IDX=/share/Passer/data/reference/house_sparrow_ref.fa.fai
# WINDOW=10000000
# OUTPUT=sparrow_genome_windows.list

# make genome size file
cut -f 1-2 ${REF_IDX} > genome_size.txt

# make the windows and cat together
bedtools makewindows -g genome_size.txt -w ${WINDOW} | \
grep -v ${SCAFFOLD_NAME} | awk '{print $1":"$2"-"$3}' \
> $OUTPUT

# GATK will fail if the coords include 0, so edit to start from 1
sed -i_bak 's/:0-/:1-/g' $OUTPUT

# creates 115 genome windows

# for scaffolds
bedtools makewindows -g genome_size.txt -w ${WINDOW} | \
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
split -d --lines=${lines_per_file} scaffolds.list2 scaffolds:

# add scafs to windows list only if any exist
n_scaffolds=$(ls scaffolds:* | wc -l)

if [ ${n_scaffolds} -gt 0 ]
then
    for i in scaffolds:*; do echo $i; done >> $OUTPUT
fi
