#!/bin/bash

# Introdução ao Shell Script no Linux
# https://www.devmedia.com.br/introducao-ao-shell-script-no-linux/25778

# echo "Seu nome de usuário é:"
# whoami
# echo "Info de hora atual e tempo que o computador está ligado:"
# uptime
# echo "O script está executando do diretório:"
# pwd

#Este é um comentário
#Este é outro comentário
#  echo "Este script contém comentários."

# site=www.devmedia.com.br
# meu_numero_favorito=13
# _cidade="Porto Alegre"
# echo "Um ótimo site para você aprender a programar e se manter atualizado é: $site"
# echo "Meu número favorito é: $meu_numero_favorito"
# echo "Minha cidade natal é: $_cidade"

# nome=fernanda
# echo "O nome da variável é \$nome"

# system_info=`df -h` # Também poderia ser system_info=$(df -h)
# echo "$system_info"

# echo "Qual o nome de uma de suas músicas favoritas?"
# read nome_musica;
# echo "Você gosta de ouvir $nome_musica"

# echo "Digite um número qualquer:"
# read numero;
# if [ "$numero" -gt 20 ];
# then
#   echo "Este número é maior que 20!"
# fi

# echo "Digite um número qualquer:"
#   read numero;
#   if [ "$numero" -ge 0 ];
#    then
#     echo "O número $numero é positivo!"
#   else
#    echo "O número $numero é negativo!"
# fi

# echo "Selecione uma opção:"
# echo "1 - Exibir data e hora do sistema"
# echo "2 - Exibir o resultado da divisão 10/2"
# echo "3 - Exibir uma mensagem"                                                        
# read opcao;                                                                         
# if [ $opcao == "1" ];
# then
#   data=$(date +"%T, %d/%m/%y, %A")
#   echo "$data"
# elif [ $opcao == "2" ];
# then
#   result=$((10/2))
#   echo "divisao de 10/2 = $result"
# elif [ $opcao == "3" ];
# then
#   echo "Informe o seu nome:"
#   read nome;
#   echo "Bem-vindo ao mundo do shell script, $nome!"
# fi

# echo "Selecione uma opção:"
#   echo "1 - Exibir data e hora do sistema"
#   echo "2 - Exibir o resultado da divisão 10/2"
#   echo "3 - Exibir uma mensagem"                                                      
#  read opcao;
#   case $opcao in
#    "1")
#       data=$(date +"%T, %d/%m/%y, %A")
#       echo "$data"
#       ;;
#    "2")
#      result=$((10/2))
#      echo "divisao de 10/2 = $result"
#    ;;
#    "3")
#     echo "Informe o seu nome:"
#     read nome;
#     echo "Bem-vindo ao mundo do shell script, $nome!"
#   ;;
# esac

# echo "Testando o loop for"
# for i in {10..0};
#   do
#     echo "$i"
# done

# echo "Testando o comando seq"
# for i in $(seq 1 100);
#   do
#     echo "$i"
# done

# echo "Testando o comando seq"
# for i in $(seq 1 5 100);
#   do
#     echo "$i"
# done

# echo "Informe o que você quiser, -1 para sair"
# read dado;
# while [ $dado != "-1" ];
#   do
#     echo "Você digitou $dado"
#     read dado;
# done

# echo "Informe até que valor positivo e maior que zero contar:"
# read valor;
# i=1
# while [ $i -le $valor ];
#   do
#     echo "$i"
#     ((i=$i+1))  
# done

# FUNÇÃO

# main()
# {
#   echo "Escolha uma opção:"
#   echo "1 - Esvaziar a lixeira"
#   echo "2 - Calcular fatorial"
#   read opcao;
#   case $opcao in
#     "1")
#       esvaziar_lixeira
#       ;;
#     "2")
#       calcular_fatorial
#       ;;
#   esac
# }

# esvaziar_lixeira()
# {
#   echo "Esvaziando a lixeira..."
#   path="${HOME}/.local/share/Trash/files"  
#   cd "$path"
#   for file in *
#     do
#       rm -rf "$file"
#   done
#   echo "Done!"
# }

# calcular_fatorial()
# {
#   echo "Informe um número:"
#   read numero;
#   i=1
#   fat=1
#   while [ $i -le $numero ]
#     do
#       fat=$(($fat*$i))
#       i=$(($i+1))
#   done
#   echo "fatorial de $numero é $fat"
# }
# main

if [ $# -lt 1 ];
  then
  echo "Precisa fornecer pelo menos 1 argumento!"
  exit 1
fi
echo "Número de argumentos passados: $#"
i=0
for argumento in $*
  do
  i=$(($i+1))
  echo "Argumento $i passado: $argumento"
done