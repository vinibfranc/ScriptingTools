from Bio import SeqIO
count = 0
for rec in SeqIO.parse("sequences/SRR020192.fastq", "fastq"):
    count += 1
print("%i reads" % count)

good_reads = (rec for rec in \
              SeqIO.parse("sequences/SRR020192.fastq", "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 20)
count = SeqIO.write(good_reads, "outputs/good_quality.fastq", "fastq")
print("Saved %i reads" % count)

'''
More information is available for trimming primer and adaptor sequences
'''