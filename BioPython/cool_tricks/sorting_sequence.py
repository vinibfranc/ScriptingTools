from Bio import SeqIO
# Get the lengths and ids, and sort on length
len_and_ids = sorted((len(rec), rec.id) for rec in
                     SeqIO.parse("ls_orchid.fasta", "fasta"))
ids = reversed([id for (length, id) in len_and_ids])
del len_and_ids  # free this memory
record_index = SeqIO.index("sequences/ls_orchid.fasta", "fasta")
records = (record_index[id] for id in ids)
SeqIO.write(records, "outputs/sorted.fasta", "fasta")