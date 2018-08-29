# ------------------------------------------------
# running the analysis for spouse gwas stuff
# ------------------------------------------------
rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

dat <- read_delim("data/differences_dat.txt", delim = "\t")

snps <- grep("rs[0-9]", colnames(dat), value = T)

# regress Height on SNPs
i=snps[1]
comp_res <- list()
sum_res <- list()

vars <- apply(dat[, snps], 2, function(x) {var(x, na.rm = T)})
vars[vars == 0] # rs9825951_A has no variance - removing

dat <- dplyr::select(dat, -one_of(names(vars[vars==0])))
snps <- grep("rs[0-9]", colnames(dat), value = T)
for (i in snps) {
	fom <- as.formula(paste0("height_diff ~ ", i, " + age_diff + Sex"))
	temp <- lm(fom, data = dat)
	x <- summarise_lm(temp, "height", i)
	comp_res[[i]] <- x
	sum_res[[i]] <- x$summary_dat
}

res <- do.call(rbind, sum_res)
res <- as.data.frame(res) %>%
	rownames_to_column(var = "snp") %>%
	mutate(snp = gsub("_[ACTG]", "", snp))

res_sig <- res %>%
	dplyr::filter(p < 0.05)

# read in actual gwas for comparison
height_gwas <- read_table("data/output/Height-Height.conventional.Height.assoc.linear")
height_gwas <- height_gwas %>%
	mutate(snp = SNP) %>%
	mutate(estimate = BETA) %>%
	mutate(se = SE) %>%
	mutate(p = P) %>%
	mutate(CI_low = gsub(" .*", "", head(height_gwas[["L95      U95"]]))) %>%
	mutate(CI_high = gsub("*. ", "", head(height_gwas[["L95      U95"]])))  

gsub("*", "", head(height_gwas[["L95      U95"]]))

# join data together - probs do it so they're side by side then remove everything but beta, se, p, l95, u95 for both 
head(fin_dat)

fin_dat <- res_sig %>%
	left_join(height_gwas, by = c("snp" = "SNP"))



# use a plot to compare - try forest?

dim(res)

x <- list(x_x = data.frame(x_x_x = 1:10))
y <- list(y_y = data.frame(y_y_y = 1:10))

z <- list(x, y)


fom <- as.formula(paste0("hdiff ~ x + DoB + Sex")
lm_res <- lm(height_diff ~ rs425277_T + age_diff + Sex, data = dat)


str(x)
lm_res <- lapply(test_diff[, -c(1:13)], function(x) {lm(fom, data = dat)})

