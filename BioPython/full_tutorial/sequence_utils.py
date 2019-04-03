from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

my_seq = Seq("GATCG", IUPAC.unambiguous_dna)
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

print("Seq length: ", len(my_seq))
print(my_seq[0]) #first letter
print(my_seq[2]) #third letter
print(my_seq[-1]) #last letter

# Non-overlapping count
print("AAAA".count("AA"))
print(Seq("AAAA").count("AA"))

# Single letter count
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
print(len(my_seq))
print(my_seq.count("G"))
# GC content
print(100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq))
# GC content builtin
print(GC(my_seq))

# Slicing a sequence

print(my_seq[0::3])
print(my_seq[1::3])
print(my_seq[2::3])
print(my_seq[::-1])

# Convert Seq object to string

str_seq = str(my_seq)
print(str_seq)

# Concatenating sequences

from Bio.Alphabet import generic_alphabet

protein_seq = Seq("EVRNAK", IUPAC.protein)
dna_seq = Seq("ACGT", IUPAC.unambiguous_dna)
protein_seq.alphabet = generic_alphabet
dna_seq.alphabet = generic_alphabet
print(protein_seq + dna_seq)

# Add multiple sequences together

from Bio.Alphabet import generic_dna

list_of_seqs = [Seq("ACGT", generic_dna), Seq("AACC", generic_dna), Seq("GGTT", generic_dna)]
print(sum(list_of_seqs, Seq("", generic_dna)))

# Upper and lower cases

print(dna_seq.upper())
print(dna_seq.lower())

print(dna_seq)

print("acg" in dna_seq)
print("ACG" in dna_seq.upper())

# Complements

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
print(my_seq)
print(my_seq.complement())
print(my_seq.reverse_complement())

# Transcription

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
print(coding_dna)
template_dna = coding_dna.reverse_complement()
print(template_dna)

messenger_rna = coding_dna.transcribe()
print(messenger_rna)
print(messenger_rna.back_transcribe())

# Translation from messenger rna

print(messenger_rna.translate())

# Translation from coding DNA

from Bio.Alphabet import IUPAC

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
print(coding_dna)
print(coding_dna.translate())

# Translation tables

from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
print(standard_table)
print(mito_table)

# Comparing seeq objects

seq1 = Seq("ACGT", IUPAC.unambiguous_dna)
seq2 = Seq("ACGT", IUPAC.ambiguous_dna)
print(str(seq1) == str(seq2))
print(str(seq1) != str(seq1))