#!/bin/bash

#bgen file
bgen_pattern=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_imp_chrCHROM_v2.bgen
#bgi file
bgen_index_pattern=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_bgi_chrCHROM_v2.bgi
#List of SNPs you want to extract
snp_list=~/spouse_gwas/data/heightSNPs.txt

# Extract using bgen to temporary files
temp_prefix=temp_genotypes
for chrom in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}; do
    inbgen=${bgen_pattern/CHROM/$chrom}
    inbgenidx=${bgen_index_pattern/CHROM/$chrom}
    ~/programs/gavinband-bgen-0b7a2803adb5/build/apps/bgenix -g $inbgen -i $inbgenidx -incl-rsids $snp_list -v11 > $temp_prefix.$chrom.bgen
done

#Use QCtools to convert to plink format
module add apps/qctool-2.0

declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")

for i in "${arr[@]}"
do
qctool -g temp_genotypes.${i}.bgen -og Temp.SNPs.${i}.bed
done


#Also will need the idlist to link the IDs to the genotype files, they are in the same order so just paste across and then remove columns in R or unix 
paste ~/spouse_gwas/idlist file.fam > linkedfile.fam

R
input <- read.table("~/spouse_gwas/linked_height_SNPs.fam")
output <- input[c(1, 1, 4:7)]

write.table(output, "~/spouse_gwas/height_SNPs.fam", sep=" ", quote=F, row.names=F, col.names=F)