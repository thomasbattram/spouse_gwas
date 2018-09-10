# ------------------------------------------------
# generating the difference files for gwas 
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "data.table")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

dat <- fread("data/merged_dat_recode.raw", header = T)
pairs <- read.table("data/spouse_pairs_height_chd.txt", header = T)
bim <- read_bim("data/binary/merged_dat")

# ------------------------------------------------
# make genotype file
# ------------------------------------------------
# If making it for a binary trait then will need to make sure the differences end up as a binary
trait <- "chd"

res <- pairs %>%
	left_join(dat, by = c("genID" = "FID"))
str(res)

snps <- grep("rs[0-9]", colnames(res), value = T)

gen <- res %>%
	arrange(Couple) %>%
	arrange(Sex) # arranging by sex ensures when looking at duplicates you pick out all males/females
head(gen[, 1:15])

if (is.binary(gen[[trait]])) {
	gen <- gen %>%
		arrange(desc(chd)) ### find a way to do this with the trait name!!!
}

dupval <- duplicated(gen$Couple) # duplicated values = couples

# Extract one from each couple
gen_c1 <- gen[dupval, ]  %>%
	arrange(Couple)
colnames(gen_c1[snps]) <- paste0(colnames(gen_c1[snps]), "_c1")
gen_c2 <- gen[!dupval, ] %>%
	arrange(Couple)
colnames(gen_c2[snps]) <- paste0(colnames(gen_c2[snps]), "_c2")
stopifnot(gen_c1$Couple == gen_c2$Couple) # Checking the couples are in the same order

gen_diff <- gen_c1
colnames(gen_diff[snps]) <- paste0(colnames(gen_diff[snps]), "_diff")

# Minus gen couple 1 from gen couple 2
gen_diff[, -c(1:13)] <- gen_c1[, -c(1:13)] - gen_c2[, -c(1:13)]
colnames(gen_diff[, -c(1:13)]) <- paste0(colnames(gen_diff[, -c(1:13)]), "_diff")
colnames(gen_diff)

summary(gen_diff[, 14:20])
# what to do with NAs???

# ------------------------------------------------
# make final difference file
# ------------------------------------------------
fin_dat <- gen_diff %>%
	mutate(trait_diff = gen_c1[[trait]] - gen_c2[[trait]]) %>%
	mutate(age_diff = gen_c1[, "DoB"] - gen_c2[, "DoB"])

write.table(fin_dat, file = paste0("data/differences_dat_", trait, ".txt"), qu = F, col = T, row = F, sep = "\t")




