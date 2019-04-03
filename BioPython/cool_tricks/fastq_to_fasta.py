from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
SeqIO.convert("outputs/good_quality.fastq", "fastq", "sequences/resulting_reads.fasta", "fasta")
SeqIO.convert("outputs/good_quality.fastq", "fastq", "sequences/resulting_reads.qual", "qual")

for record in PairedFastaQualIterator(open("sequences/resulting_reads.fasta"), open("sequences/resulting_reads.qual")):
   print(record)