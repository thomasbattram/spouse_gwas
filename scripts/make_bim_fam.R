# ------------------------------------------------
# Making a fam and bim file for the spouse pairs in biobank
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

pairs <- read_delim("data/spouse_pairs.txt", delim = " ")
snps <- read_delim("data/height_mr_dataset.txt", delim = " ")
# ------------------------------------------------
# Make fam file
# ------------------------------------------------

# FAM format
# FID IID PID MID SEX PHENO
# --  --  --  --  --  ---
# --  --  --  --  --  ---
# --  --  --  --  --  ---
# --  --  --  --  --  ---

fam <- pairs %>%
	mutate(FID = Couple) %>%
	mutate(IID = genID) %>%
	mutate(PID = 0) %>%
	mutate(MID = 0) %>%
	mutate(SEX = as.factor(Sex)) %>%
	mutate(PHENO = 0) %>%
	dplyr::select(FID, IID, PID, MID, SEX, PHENO)


# ------------------------------------------------
# Make bim file
# ------------------------------------------------

# BIM format
# CHR SNP/CPG GD POS A1 A2
# --	--    -- --  -- --
# --	--    -- --  -- --
# --	--    -- --  -- --
# --	--    -- --  -- --

filepath <- "/panfs/panasas01/shared/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/filtered/bestguess/maf0.01_info0.8/combined"

als_bim <- read_delim(paste0(filepath, "/data.bim"), delim = "\t", col_names = F)
colnames(als_bim) <- c("CHR", "SNP", "GD", "POS", "A1", "A2")
bim <- snps %>%
	left_join(als_bim) %>% 
	mutate(A1 = effect_allele.exposure) %>%
	mutate(A2 = other_allele.exposure) %>%
	dplyr::select(CHR, SNP, GD, POS, A1, A2)

bim

for (i in 1:22) {
	
}


