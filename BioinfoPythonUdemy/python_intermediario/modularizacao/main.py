import aleatorio as a
import media as m

lista = a.geraListaInteiros(10)
print(lista)

media = m.media(lista)
mediana = m.mediana(lista)
moda = m.moda(lista)

print("Média da lista: "+str(media))
print("Mediana da lista: "+str(mediana))
print("Moda da lista: "+str(moda))