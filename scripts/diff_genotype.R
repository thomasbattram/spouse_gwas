# ------------------------------------------------
# generating the genotype for gwas 
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

dat <- read.table("data/merged_dat_recode.raw", header = T)
pairs <- read.table("data/spouse_pairs.txt", header = T)

# ------------------------------------------------
# make genotype file
# ------------------------------------------------

res <- pairs %>%
	left_join(dat, by = c("genID" = "FID"))
str(res)

snps <- grep("rs[0-9]", colnames(res), value = T)

test <- res[1:20, ] %>%
	arrange(Sex) # arranging by sex ensures when looking at duplicates you pick out all males/females

dupval <- duplicated(test$Couple) # duplicated values = couples

# Extract one from each couple
test_c1 <- test[dupval, ] 
colnames(test_c1[snps]) <- paste0(colnames(test_c1[snps]), "_c1")
test_c2 <- test[!dupval, ]
colnames(test_c2[snps]) <- paste0(colnames(test_c2[snps]), "_c2")
stopifnot(test_c1$Couple == test_c2$Couple)

test_diff <- test[dupval, ]
colnames(test_diff[snps]) <- paste0(colnames(test_diff[snps]), "_diff")

# Minus test couple 1 from test couple 2
test_diff[, -c(1:11)] <- test_c1[,-c(1:11)] - test_c2[,-c(1:11)]
colnames(test_diff[, -c(1:11)]) <- paste0(colnames(test_diff[, -c(1:11)]), "_diff")
colnames(test_diff)
# START FROM HERE!!!!!!
summary(test_diff[,12:20])
test_diff[, -c(1:11)] <- 


fin_dat <- list()
i=1
for (i in 1:length(snps)) {
	exlude_snp <- snps[-i]
	temp <- test %>%
		dplyr::select(-one_of(exlude_snp)) %>%
		spread(Couple)
}

test <- res %>%
	gather(variable, value, -(phenID:PHENOTYPE)) %>%
	unite(temp, Couple, variable) %>%
	spread(temp, value)

test[1:10, 1:20]
dim(test)
head(test)


df <- data.frame(month=rep(1:3,2),
                 student=rep(c("Amy", "Bob"), each=3),
                 A=c(9, 7, 6, 8, 6, 9),
                 B=c(6, 7, 8, 5, 6, 7))


df %>% 
  gather(variable, value, -(month:student)) %>%
  unite(temp, student, variable) %>%
  spread(temp, value)
