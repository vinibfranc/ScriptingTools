from Bio import SeqIO

'''Reading sequence files'''

# FASTA

# for seq_record in SeqIO.parse("sequences/ls_orchid.fasta", "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

# GenBank

# for seq_record in SeqIO.parse("sequences/ls_orchid.gbk", "genbank"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

'''Iterating over the records in a sequence file'''

# record_iterator = SeqIO.parse("sequences/ls_orchid.fasta", "fasta")

# first_record = next(record_iterator)
# print(first_record.id)
# print(first_record.description)

# second_record = next(record_iterator)
# print(second_record.id)
# print(second_record.description)

'''Getting a list of the records in a sequence file'''

# records = list(SeqIO.parse("sequences/ls_orchid.gbk", "genbank"))

# print("Found %i records" % len(records))

# print("The last record")
# last_record = records[-1] #using Python's list tricks
# print(last_record.id)
# print(repr(last_record.seq))
# print(len(last_record))

# print("The first record")
# first_record = records[0] #remember, Python counts from zero
# print(first_record.id)
# print(repr(first_record.seq))
# print(len(first_record))

'''Extracting data'''

# record_iterator = SeqIO.parse("sequences/ls_orchid.gbk", "genbank")
# first_record = next(record_iterator)
# print(first_record)

'''Organism information'''

# all_species = []
# for seq_record in SeqIO.parse("sequences/ls_orchid.gbk", "genbank"):
#     all_species.append(seq_record.annotations["organism"])
# print(all_species)

# all_species = []
# for seq_record in SeqIO.parse("sequences/ls_orchid.fasta", "fasta"):
#     all_species.append(seq_record.description.split()[1])
# print(all_species)

'''Parsing sequences from compressed files'''

# import gzip
# with gzip.open("sequences/ls_orchid.gbk.gz", "rt") as handle:
#     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))

'''Parsing GenBank records from the net'''

# from Bio import Entrez
# Entrez.email = "viniciusfr@ufcspa.edu.br"
# with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="6273291") as handle:
#     seq_record = SeqIO.read(handle, "fasta")
# print("%s with %i features" % (seq_record.id, len(seq_record.features)))

'''Sequence files as Dictionaries – In memory'''

# orchid_dict = SeqIO.to_dict(SeqIO.parse("sequences/ls_orchid.gbk", "genbank"))
# print(len(orchid_dict))
# print(list(orchid_dict.keys()))
# print(list(orchid_dict.values()))

# seq_record = orchid_dict["Z78475.1"]
# print(seq_record.description)
# print(repr(seq_record.seq))

'''Sequence files as Dictionaries – Database indexed files'''

# Run it on terminal

# rsync -avP "ftp.ncbi.nih.gov::genbank/gbvrl*.seq.gz" .
# gunzip gbvrl*.seq.gz
# curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl1.seq.gz
# curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl2.seq.gz
# curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl3.seq.gz
# curl -O ftp://ftp.ncbi.nih.gov/genbank/gbvrl4.seq.gz
# gunzip gbvrl*.seq.gz

# Extracting info

# import glob
# files = glob.glob("gbvrl*.seq")
# print("%i files to index" % len(files))
# gb_vrl = SeqIO.index_db("gbvrl.idx", files, "genbank")
# print("%i sequences indexed" % len(gb_vrl))

# Raw data

# print(gb_vrl.get_raw("AB811634.1"))

'''Low level FASTA and FASTQ parsers'''

# FASTA

# from Bio.SeqIO.FastaIO import SimpleFastaParser
# count = 0
# total_len = 0
# with open("sequences/ls_orchid.fasta") as in_handle:
#     for title, seq in SimpleFastaParser(in_handle):
#         count += 1
#         total_len += len(seq)
# print("%i records with total sequence length %i" % (count, total_len))

# FASTQ

# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# count = 0
# total_len = 0
# with open("example.fastq") as in_handle:
#     for title, seq, qual in FastqGeneralIterator(in_handle):
#         count += 1
#         total_len += len(seq)
# print("%i records with total sequence length %i" % (count, total_len))




