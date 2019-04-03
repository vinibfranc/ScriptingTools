'''
The SeqRecord class itself is quite simple, and offers the following information as attributes:

.seq
– The sequence itself, typically a Seq object.
.id
– The primary ID used to identify the sequence – a string. In most cases this is something like an accession number.
.name
– A “common” name/id for the sequence – a string. In some cases this will be the same as the accession number, but it could also be a clone name. I think of this as being analogous to the LOCUS id in a GenBank record.
.description
– A human readable description or expressive name for the sequence – a string.
.letter_annotations
– Holds per-letter-annotations using a (restricted) dictionary of additional information about the letters in the sequence. The keys are the name of the information, and the information is contained in the value as a Python sequence (i.e. a list, tuple or string) with the same length as the sequence itself. This is often used for quality scores (e.g. Section 20.1.6) or secondary structure information (e.g. from Stockholm/PFAM alignment files).
.annotations
– A dictionary of additional information about the sequence. The keys are the name of the information, and the information is contained in the value. This allows the addition of more “unstructured” information to the sequence.
.features
– A list of SeqFeature objects with more structured information about the features on a sequence (e.g. position of genes on a genome, or domains on a protein sequence). The structure of sequence features is described below in Section 4.3.
.dbxrefs
- A list of database cross-references as strings.

'''

# Usual Record
from Bio.Seq import Seq
simple_seq = Seq("GATC")

# Creating SeqRecord

from Bio.SeqRecord import SeqRecord
simple_seq_r = SeqRecord(simple_seq)
print(simple_seq_r.id)
simple_seq_r.id = "AC12345"
simple_seq_r.description = "Made up sequence I wish I could write a paper about"
print(simple_seq_r.description)
print(simple_seq_r.seq)
simple_seq_r = SeqRecord(simple_seq, id="AC12345")
print(simple_seq_r)
simple_seq_r.annotations["evidence"] = "None. I just made it up."
print(simple_seq_r.annotations)
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
print(simple_seq_r.letter_annotations)
print(simple_seq_r.letter_annotations["phred_quality"])

# SeqRecord from FASTA
# from Bio import SeqIO
# record = SeqIO.read("NC_005816.fna", "fasta")
# print(record)
# print(record.seq)
# print(record.id)
# print(record.name)
# print(record.description)
# print(record.dbxrefs)
# print(record.annotations)
# print(record.letter_annotations)
# print(record.features)

# Comparison

record1 = SeqRecord(Seq("ACGT"), id="test")
record2 = SeqRecord(Seq("ACGT"), id="test")
# print(record1 == record2)
print(record1.id == record2.id)


