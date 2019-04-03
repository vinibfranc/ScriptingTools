from Bio import SeqIO

fq_dict = SeqIO.index("sequences/SRR020192.fastq", "fastq")
print(len(fq_dict))
print(fq_dict.keys())
print(fq_dict["SRR020192.23186"].seq)