# ------------------------------------------------
# Making a fam file 
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

# FAM format
# FID IID PID MID SEX PHENO
# --  --  --  --  --  ---
# --  --  --  --  --  ---
# --  --  --  --  --  ---
# --  --  --  --  --  ---

fam <- read_delim("data/binary/merged_dat.fam", delim = " ", col_names = F)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

linker <- read_delim("data/idlist", delim = " ", col_names = F)

fam[,c(1,2)] <- linker[[1]]

write.table(fam, "data/binary/merged_dat.fam", col.names = F, row.names = F, quote = F)


