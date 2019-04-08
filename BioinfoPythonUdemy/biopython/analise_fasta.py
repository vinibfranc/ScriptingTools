from Bio import SeqIO

for fasta in SeqIO.parse("exemplo.fasta", "fasta"):
    print(fasta.id)
    print(fasta.seq)
