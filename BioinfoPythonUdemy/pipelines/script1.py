# OS

import os

os.system("python3 script2.py")

# System não permite que haja interação entre resultado dos scripts
# Solução: redirecionar para arquivo

variavel = os.system("python3 script2.py > output.txt")

arquivo = open("output.txt").read()

print("A mensagem executada por script2.py foi: "+arquivo)

# subprocess

import subprocess

print("PARTE 1")
variavel = subprocess.call(["python", "script2.py"])

print("PARTE 2")
print(variavel)

print("PARTE 1")
variavel = subprocess.check_output(["python", "script2.py"])

print("PARTE 2")
print(variavel)

