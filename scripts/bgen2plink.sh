#!/bin/bash

cd ~/spouse_gwas/

# get the snps from the mr base output
awk 'NR!=1{print $1}' data/height_mr_dataset.txt > data/height_SNPs.txt

#bgen file
bgen_pattern=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_imp_chrCHROM_v2.bgen
#bgi file
bgen_index_pattern=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_bgi_chrCHROM_v2.bgi
#List of SNPs you want to extract
snp_list=~/spouse_gwas/data/height_SNPs.txt

# Extract using bgen to temporary files
temp_prefix=temp_genotypes
for chrom in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}; do
    inbgen=${bgen_pattern/CHROM/$chrom}
    inbgenidx=${bgen_index_pattern/CHROM/$chrom}
    ~/programs/gavinband-bgen-0b7a2803adb5/build/apps/bgenix -g $inbgen -i $inbgenidx -incl-rsids $snp_list -v11 > data/bgen/$temp_prefix.$chrom.bgen
done

#Use QCtools to convert to plink format
module rm apps/qctool-1.4
module add apps/qctool-2.0

declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")

for i in "${arr[@]}"
do
qctool -g data/bgen/temp_genotypes.${i}.bgen -og data/binary/temp_SNPs.${i}.bed
done

gd=/panfs/panasas01/shared-biobank/data/bestguess

grep -hf $snp_list $gd/*.bim > data/snp_pos.txt
wc -l data/snp_pos.txt # many more than 382... - check in R if this can be sorted


#Also will need the idlist to link the IDs to the genotype files, they are in the same order so just paste across and then remove columns in R or unix 
Rscript scripts/make_fam.R
