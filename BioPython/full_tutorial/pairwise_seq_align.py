'''Create the alignment'''

from Bio import Align
aligner = Align.PairwiseAligner()

'''Calculating score'''

# Default = global alignment

# seq1 = "GAACT"
# seq2 = "GAT"
# score = aligner.score(seq1, seq2)
# print(score)

# alignments = aligner.align(seq1, seq2)
# for alignment in alignments:
#     print(alignment)

'''Local alignment'''

# aligner.mode = 'local'
# seq1 = "AGAACTC"
# seq2 = "GAACT"
# score = aligner.score(seq1, seq2)
# print(score)

# alignments = aligner.align(seq1, seq2)
# for alignment in alignments:
#     print(alignment)

# print(aligner)
# print(aligner.algorithm)
# print(aligner.epsilon)

'''Match and mismatch scores'''

# aligner = Align.PairwiseAligner()
# print(aligner.match)
# print(aligner.mismatch)
# score = aligner.score("AAA","AAC")
# print(score)
# aligner.match = 2.0
# score = aligner.score("AAA","AAC")
# print(score)

# Specifying a substitution matrix

# matrix = {('A','A'): 1.0, ('A','B'): -1.0, ('B','B'): 2.0}
# aligner.substitution_matrix = matrix
# aligner.substitution_matrix[('A','A')] = 5.0 # does nothing
# print(aligner.substitution_matrix[('A','A')])

'''More details in: 
    General gap scores, 
    Iterating over alignments,
    Alignment objects
'''