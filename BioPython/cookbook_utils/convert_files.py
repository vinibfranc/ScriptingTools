from Bio import SeqIO

with open("sequences/cor6_6.gb", "rU") as input_handle:
    with open("results/cor6_6.fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")

print("Converted %i records" % count)