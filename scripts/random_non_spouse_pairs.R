# ------------------------------------------------
# Comparison of within spouse gwas and random non-spouse pairs gwas
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "data.table", "ggExtra")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

dat <- fread("data/merged_dat_recode.raw", header = T)
pairs <- read.table("data/spouse_pairs_height_chd.txt", header = T)
pairs$height <- pairs$Height
bim <- read_bim("data/binary/merged_dat")

args <- commandArgs(trailingOnly = T)
num <- args[1]
traits <- args[2]
outcome <- args[3]
set.seed(num)
traits <- unlist(strsplit(traits, ","))
print(traits)
# ------------------------------------------------
# make genotype file
# ------------------------------------------------
# If making it for a binary trait then will need to make sure the differences end up as a binary
# traits <- c("chd", "height")
# outcome <- "chd"

res <- pairs %>%
	left_join(dat, by = c("genID" = "FID"))
str(res)

snps <- grep("rs[0-9]", colnames(res), value = T)

gen <- res %>%
	arrange(Couple) %>%
	arrange(Sex) # arranging by sex ensures when looking at duplicates you pick out all males/females
head(gen[, 1:15])
if (is.binary(gen[[outcome]])) {
	temp <- quos(!! sym(outcome))
	gen <- gen %>%
		arrange(desc(!!! temp)) 
}

gen_random_pairs <- function(dat, num) {
	x <- sample(num)
	y <- sample(num)
	temp <- dat %>%
		mutate(random_pair = NA)

	temp[temp$Sex == 1, "random_pair"] <- x
	temp[temp$Sex == 0, "random_pair"] <- y
	return(temp)
}

gen_dat <- gen_random_pairs(gen, nrow(gen)/2)
couple_tab <- gen_dat %>%
	dplyr::select(Couple, Sex, random_pair) %>%
	spread(key = Sex, value = random_pair)
same_pair <- couple_tab[couple_tab[["0"]] == couple_tab[["1"]], "1"]
print(paste0("number of pairs flipping = ", length(same_pair)))
while (length(same_pair) > 0) {
	rand_swap <- sample(gen_dat$random_pair, 1)
	pos1 <- which(gen_dat$random_pair == rand_swap)[1]
	pos2 <- which(gen_dat$random_pair == same_pair[1])[1]
	gen_dat$random_pair <- replace(gen_dat$random_pair, c(pos1, pos2), gen_dat[c(pos2, pos1), "random_pair"])
	couple_tab <- gen_dat %>%
		dplyr::select(Couple, Sex, random_pair) %>%
		spread(key = Sex, value = random_pair)
	same_pair <- couple_tab[couple_tab[["0"]] == couple_tab[["1"]], "1"]
}	
stopifnot(max(table(gen_dat$random_pair)) == 2)

dupval <- duplicated(gen_dat$random_pair)  # duplicated values = pairs
gen_c1 <- gen_dat[dupval, ]  %>%
	arrange(random_pair)
colnames(gen_c1[snps]) <- paste0(colnames(gen_c1[snps]), "_c1")
gen_c2 <- gen_dat[!dupval, ]  %>%
	arrange(random_pair)
colnames(gen_c2[snps]) <- paste0(colnames(gen_c2[snps]), "_c2")
stopifnot(gen_c1$random_pair == gen_c2$random_pair) # Checking the couples are in the same order

gen_diff <- gen_c1
colnames(gen_diff[snps]) <- paste0(colnames(gen_diff[snps]), "_diff")

# Minus gen couple 1 from gen couple 2
gen_diff[, -c(1:13)] <- gen_c1[, -c(1:13)] - gen_c2[, -c(1:13)]
colnames(gen_diff[, -c(1:13)]) <- paste0(colnames(gen_diff[, -c(1:13)]), "_diff")
colnames(gen_diff)

summary(gen_diff[, 14:20])

# remove SNPs that are missing in over 10% of individuals
na_num <- apply(gen_diff, 2, function(x) {sum(is.na(x))})
summary(na_num)
rm_dat <- names(na_num[na_num > nrow(gen_diff) / 10])
length(rm_dat)
gen_diff <- gen_diff %>%
	dplyr::select(-one_of(rm_dat))
if (any(traits %in% rm_dat)) {
	stop("TOO MUCH MISSING DATA IN TRAIT(S)")
}

# remove individuals that have missing data for over 10% of snps
na_num <- apply(gen_diff, 1, function(x) {sum(is.na(x))})
summary(na_num)
names(na_num) <- 1:nrow(gen_diff)
snp_start <- min(grep("rs[0-9]", colnames(gen_diff)))
rm_dat <- as.numeric(names(na_num[na_num > (ncol(gen_diff) - snp_start) / 10]))
length(rm_dat)
if (length(rm_dat) != 0) {
	gen_diff <- gen_diff[-rm_dat, ]
}
print("starting final bit")
# ------------------------------------------------
# make final difference file
# ------------------------------------------------
fin_dat <- gen_diff %>%
	# mutate(trait_diff = gen_c1[[trait]] - gen_c2[[trait]]) %>%
	mutate(age_diff = gen_c1[, "DoB"] - gen_c2[, "DoB"]) %>%
	mutate(sex_diff = gen_c1[, "Sex"] - gen_c2[, "Sex"])

for (i in traits) {
	var_nam <- paste0(i, "_diff")
	fin_dat[[var_nam]] <- gen_c1[[i]] - gen_c2[[i]]
}

write.table(fin_dat, file = paste0("data/random_pairs/pair", num, "_differences_dat_", paste(traits, collapse = "_"), ".txt"), qu = F, col = T, row = F, sep = "\t")
