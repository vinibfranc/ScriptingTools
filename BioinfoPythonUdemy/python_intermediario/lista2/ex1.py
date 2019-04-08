import re

sequencia1 = input("Digite uma sequência: ")
sequencia2 = input("Digite outra sequência: ")

# if sequencia1 == sequencia2:
#     print("Sequências iguais!")
# else:
#     print("Sequências diferentes!")

busca = re.match(sequencia1, sequencia2)

if busca:
    print("Sequências idênticas")
    print(busca.group())
else:
    print("Sequências diferentes")