import sys
from Bio import SeqIO
SeqIO.convert(sys.stdin, "fastq-solexa", sys.stdout, "fastq")

# Run on terminal: 
# gunzip -c some_solexa.fastq.gz | python solexa2sanger_fq.py

# Redirection: python solexa2sanger_fq.py < some_solexa.fastq > some_phred.fastq