from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# DNA Sequence

my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
print(my_seq)
print(my_seq.alphabet)

# Protein Sequence
my_seq = Seq("AGTACACTGGT", IUPAC.protein)
print(my_seq)
print(my_seq.alphabet)