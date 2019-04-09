import os

os.system("blastp -query a.fasta -subject b.fasta -outfmt 6 > resultado.txt")

linhas = open("resultado.txt").readlines()

for linha in linhas:
    print(linha)