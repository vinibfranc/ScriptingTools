nome = input("Digite o nome do arquivo que você deseja abrir: ")

arquivo = open(nome)
linhas = arquivo.readlines()

for linha in linhas:
    print(linha.strip())

