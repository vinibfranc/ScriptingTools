from Bio.Seq import Seq

minha_seq = Seq("ATCGATGCA")
print(minha_seq)

# seq complementar
seq_complementar = minha_seq.complement()
print(seq_complementar)

# seq reversa complementar

seq_complementar_reversa = minha_seq.reverse_complement()
print(seq_complementar_reversa)

# Transcrição

seq_rna = minha_seq.transcribe()
print(seq_rna)

# Retornar para DNA

seq_dna = seq_rna.back_transcribe()
print(seq_dna)

# Tradução

seq_proteina = seq_rna.translate()
print(seq_proteina)

seq_protenina_dna = seq_dna.translate()
print(seq_protenina_dna)