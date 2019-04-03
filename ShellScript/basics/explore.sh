#!/bin/bash
cat my_file.txt
echo "------------"
head my_file.txt
echo "------------"
head -n 5 my_file.txt
echo "------------"
less -S my_file.txt
echo "------------"
cat my_file.txt | cut -f 2
echo "------------"
cat my_file.txt | awk '{$1=$1;print}' OFS='\t' | cut -f 2
echo "------------"
cat my_file.txt | cut -f 2 | sort
echo "------------"
cat my_file.txt | cut -f 2 | uniq -c 
echo "------------"
cat my_file.txt | cut -f 2 | sort | uniq -c | sort -k1,1nr
echo "------------"