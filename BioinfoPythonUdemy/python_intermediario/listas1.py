minha_lista = ["abacaxi", "melancia", "abacate"]
minha_lista2 = [1, 2, 3, 4, 5]
minha_lista3 = ["abacaxi", 2, 9.89, True]

for item in minha_lista:
    print(item)

minha_lista.append("limao")
print(minha_lista)

# Verificar se elemento pertence a lista

if 3 in minha_lista2:
    print("3 está na lista2")

del minha_lista[2:]
print(minha_lista)

minha_lista4 = []

minha_lista4.append(57)
print(minha_lista4)
