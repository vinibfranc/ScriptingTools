from Bio.PDB import *

# SMCRA

parser = PDBParser()
estrutura = parser.get_structure("1BGA", "1bga.pdb")

# modelo = estrutura[0]
# cadeia_A = modelo['A']
# residuo_100 = cadeia_A[100]
# ca_100 = residuo_100['CA']

atomo_100 = estrutura[0]['A'][100]['CA']
atomo_101 = estrutura[0]['A'][101]['CA']

# distância euclidiana entre dois átomos

distancia = atomo_101 - atomo_100
print(distancia)
