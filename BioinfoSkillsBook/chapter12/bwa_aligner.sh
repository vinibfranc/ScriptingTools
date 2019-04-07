#!/bin/bash
set -e
set -u
set -o pipefail
# specify the input samples file, where the third
# column is the path to each sample FASTQ file
sample_info=samples.txt
# our reference
reference=zmays_AGPv3.20.fa
# create a Bash array from the first column, which are
# sample names. Because there are duplicate sample names
# (one for each read pair), we call uniq
sample_names=($(cut -f 1 "$sample_info" | uniq))
for sample in ${sample_names[@]}
do
# create an output file from the sample name
    results_file="${sample}.sam"
    bwa mem $reference ${sample}_R1.fastq ${sample}_R2.fastq \
> $results_file
done