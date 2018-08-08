# ------------------------------------------------
# Making a fam and bim file for the spouse pairs in biobank
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

pairs <- read_delim("data/spouse_pairs.txt", delim = " ")
snps <- read_delim("data/height_mr_dataset.txt", delim = " ")
snps2 <- read_delim("data/snp_pos.txt", delim = "\t", col_names = F)
colnames(snps2) <- c("CHR", "SNP", "GD", "POS", "A1", "A2")

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

bim <- snps %>%
	left_join(snps2) %>%
	mutate(A1 = effect_allele.exposure) %>%
	mutate(A2 = other_allele.exposure) %>%
	dplyr::select(CHR, SNP, GD, POS, A1, A2)

for (i in 1:22) {
	temp_bim <- bim %>%
		dplyr::filter(CHR == i)

	write.table(temp_bim, file = paste0("data/binary/temp_SNPs.", i, ".bim"), quote = F, col.names = F, row.names = F)
	write.table(fam, file = paste0("data/binary/temp_SNPs.", i, ".fam"), quote = F, col.names = F, row.names = F)
}

