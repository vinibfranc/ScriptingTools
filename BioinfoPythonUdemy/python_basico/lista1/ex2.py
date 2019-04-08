nota1 = float(input("Digite a primeira nota: "))
nota2= float(input("Digite a segunda nota: "))
media = (nota1 + nota2) / 2
print("Sua média é: ", media)

if(media >= 0 and media < 6):
    print("Reprovado")
elif(media >= 6 and media <= 10):
    print("Aprovado")
else:
    print("Média inválida")