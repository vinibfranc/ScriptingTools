arquivo = open("meu_arquivo.txt")

linhas = arquivo.readlines()

print(linhas)

for linha in linhas:
    print(linha)


texto = arquivo.read()

print(texto)

w = open("meu_arquivo2.txt", "w")
w.write("Esse é o meu arquivo. \nIsso fica em nova linha \n")
w.close()

w = open("meu_arquivo2.txt", "a")
w.write("Esse é o meu arquivo. \nIsso fica em nova linha")
w.close()
