idade = int(input("Digite sua idade: "))

if(idade >= 0 and idade < 18):
    print("Menor de idade!")
elif(idade >= 18 and idade <= 150):
    print("Maior de idade!")
else:
    print("Idade inválida!")
    
    