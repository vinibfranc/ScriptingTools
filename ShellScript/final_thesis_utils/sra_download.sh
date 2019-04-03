#!/bin/bash

# https://github.com/Gian77/Bash-For-Bioinformatics

# first download the accessions list of the SRA archive 
# and save them into a file.txt then run the script below

while read line; do
echo $line
fastq-dump -I --split-files $line
done < SraAccList.txt
