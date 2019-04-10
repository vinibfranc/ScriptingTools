# ****************************************************************************************
# importações
# ****************************************************************************************
import os
import requests
from Bio.PDB import *
 
 
# ****************************************************************************************
# Variáveis importantes
# ****************************************************************************************
 
sequencia = "beta_glicosidase.fasta"
diretorio = "/home/vinibfranc/Documentos/Scripting/BioinfoPythonUdemy/pipelines/lista3"
total_modelos = 2
 
 
# ****************************************************************************************
# BLAST
# ****************************************************************************************
os.system("blastp -query beta_glicosidase.fasta -subject fasta.txt -outfmt=6 > resultado.txt")
 
# Lendo o resultado do BLAST
resultado = open("resultado.txt").readlines()
maior = 0
identidade = ''
id_maior = ''
for linha in resultado:
    colunas = linha.split("\t")
    atual = float(colunas[11])
    identidade = float(colunas[2])
    if(atual > maior):
        if(identidade > 25):
            maior = atual
            id_maior = colunas[1]
            cadeia = id_maior[5] 
            id_maior = id_maior[0:4]
print("O PDB template é: "+id_maior+" (IDENTIDADE = "+str(identidade)+"% - SCORE = "+str(maior)+")")
 
# download do arquivo pdb
url = "https://files.rcsb.org/download/"+id_maior+".pdb"
r = requests.get(url)
w = open(id_maior+".pdb","w")
w.write(r.text)
w.close()
 
 
# ****************************************************************************************
# Alinhamento
# ****************************************************************************************
 
# seq = beta_glicosidase.fasta
# template = 3UR7.fasta
w = open(id_maior+".fasta","w")
w.write(">"+id_maior+"\n")
comeco = None
fim = 0
pdb = open(id_maior+".pdb").readlines()
for linha in pdb:
    if(linha[0:4] == "ATOM" and linha[21] == cadeia and linha[13:15] == 'CA'):
        resname3 = linha[17:20]
        if(comeco == None):
            comeco = int(linha[22:26])
        if(int(linha[22:26]) > fim):
            fim = int(linha[22:26])
        resname1 = Polypeptide.three_to_one(resname3)
        w.write(resname1)
w.write("\n")
w.close()
os.system("cat "+id_maior+".fasta > alinha.fasta")
os.system("cat "+sequencia+" >> alinha.fasta")
 
# run clustal-w
os.system("clustalw -infile='alinha.fasta' -output='pir'")
aln = open("alinha.pir").readlines()
new_aln = open("new_alinha.pir","w")
tipo = 0 #0 = PDB; 1 = SEQ
seq = open(sequencia)
seq_final = ""
for linha in seq:
    if(linha[0] != ">"):
        seq_final += linha.strip()
tamanho_seq = len(seq_final)
print("tamanho da seq = "+str(tamanho_seq))
for linha in aln:
    if(linha[0] == ">"):
        if(tipo == 0 and linha != "\n"):
            new_aln.write(">P1;"+id_maior+"\n")
            new_aln.write("structure:"+id_maior+":"+str(comeco)+":"+cadeia+":"+str(fim)+":"+cadeia+"::::\n")
            tipo = tipo+1
        elif(tipo == 1):
            new_aln.write(">P1;"+sequencia+"\n")
            new_aln.write("sequence:"+sequencia+":"+str(1)+":A:"+str(tamanho_seq)+":A::::")
    else:
        new_aln.write(linha)
new_aln.close()
 
 
# ****************************************************************************************
# modeller
# ****************************************************************************************
 
w = open("run.py","w")
script = "\
#-*- coding:utf-8 -*-\n\
from modeller import *\n\
from modeller.automodel import *\n\
\n\
log.verbose()\n\
env = environ()\n\
\n\
env.io.atom_files_directory = ['"+diretorio+"']\n\
\n\
env.io.hetatm = True\n\
env.io.water = True\n\
\n\
a = automodel(env, alnfile = 'new_alinha.pir', knowns = '"+id_maior+"', sequence = '"+sequencia+"')\n\
a.starting_model = 1\n\
a.ending_model = "+str(total_modelos)+"\n\
a.make() \n\
"
w.write(script)
w.close()
os.system("mod9.21 run.py")

# melhor modelo
menor = 9999999999999
menor_id = ''

for _, _, arquivos in os.walk('/home/vinibfranc/Documentos/Scripting/BioinfoPythonUdemy/pipelines/lista3'):
    
    for arquivo in arquivos:
        if(arquivo[-4:] == ".pdb"):
            try:
                linha = open(arquivo).readlines()

                score = linha[1]
                score = float(score[40:])

                if score < menor:
                    menor = score
                    menor_id = arquivo
            except:
                print("Modelo inválido")
print("O melhor modelo é: "+menor_id)

