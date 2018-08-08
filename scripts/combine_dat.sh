#!/bin/bash

cd ~/spouse_gwas/data/binary

ls > files.txt
nano files.txt # Remove files.txt from files.txt

awk '{gsub(".bim", "");print}' files.txt > file2.txt
awk '{gsub(".fam", "");print}' file2.txt > file3.txt
awk '{gsub(".bed", "");print}' file3.txt > file4.txt

plink --bfile temp_SNPs.1 --merge-list file4.txt --make-bed --out merged_dat

rm file*