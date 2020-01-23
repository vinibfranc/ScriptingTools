def dobro(x):
    return x*2

valores = [1, 2, 3, 4, 5]

lista_dobro = map(dobro, valores)
lista_dobro = list(lista_dobro)
print(lista_dobro)