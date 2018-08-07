#!/bin/bash

#bgen file
bgen_pattern=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_imp_chrCHROM_v2.bgen
#bgi file
bgen_index_pattern=/panfs/panasas01/dedicated-mrcieu/research/data/ukbiobank/_latest/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/raw_downloaded/dosage_bgen/ukb_bgi_chrCHROM_v2.bgi
#List of SNPs you want to extract
snp_list=/panfs/panasas01/sscm/lh14833/Assortative/heightSNPs.txt

# Extract using bgen to temporary files
temp_prefix=temp_genotypes
for chrom in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}; do
    inbgen=${bgen_pattern/CHROM/$chrom}
    inbgenidx=${bgen_index_pattern/CHROM/$chrom}
    /panfs/panasas01/sscm/lh14833/gavinband-bgen-798eca81c0fa/build/apps/bgenix -g $inbgen -i $inbgenidx -incl-rsids $snp_list -v11  > $temp_prefix.$chrom.bgen
done

#Use QCtools to convert to plink format
module add apps/qctool-2.0

declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")

for i in "${arr[@]}"
do
qctool -g temp_genotypes.${i}.bgen -og Temp.SNPs.${i}.bed
done


#Also will need the idlist to link the IDs to the genotype files, they are in the same order so just paste across and then remove columns in R or unix 
paste ~/Assortative/idlist file.fam > linkedfile.fam

R
input<-read.table("~/Assortative/linkedHeight.SNPs.fam")
output<-input[c(1,1,4:7)]

write.table(output, "~/Assortative/Height.SNPs.fam", sep=" ", quote=F, row.names=F, col.names=F)