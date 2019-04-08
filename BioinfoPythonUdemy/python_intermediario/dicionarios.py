meu_dicionario = {"A": "Adenina", "T": "Timina", "C": "Citosina", "G": "Guanina"}

print(meu_dicionario["A"])
print(meu_dicionario["T"])
print(meu_dicionario["C"])
print(meu_dicionario["G"])

print(meu_dicionario)

for chave in meu_dicionario:
    print(chave+" : "+meu_dicionario[chave])

for i in meu_dicionario.items():
    print(i)

for i in meu_dicionario.keys():
    print(i)

for i in meu_dicionario.values():
    print(i)