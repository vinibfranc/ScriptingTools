seq = input("Digite uma sequência: ")

arquivo = open("seq2.fasta", "w")

arquivo.write(">seq\n")
arquivo.write(seq)
arquivo.close()