# ------------------------------------------------
# running the analysis for spouse gwas stuff
# ------------------------------------------------
rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "ggExtra")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")
traits <-  c("chd","height")
dat <- read_delim(paste0("data/differences_dat_", paste(traits, collapse = "_"), ".txt"), delim = "\t")

snps <- grep("rs[0-9]", colnames(dat), value = T)

# regress height on SNPs
i=snps[1]
comp_res <- list()
sum_res <- list()

vars <- apply(dat[, snps], 2, function(x) {var(x, na.rm = T)})
vars[vars == 0] # rs9825951_A has no variance - removing

dat <- dplyr::select(dat, -one_of(names(vars[vars==0])))
snps <- grep("rs[0-9]", colnames(dat), value = T)
i=snps[1]
j=traits[1]
comp_res2 <- list()
sum_res2 <- list()
for (j in traits) {
	for (i in snps) {
		print(j)
		fom <- as.formula(paste0(j, "_diff ~ ", i, " + age_diff + Sex"))
		if (is.binary(dat[[j]])) {
			temp <- glm(fom, family = "binomial", data = dat)
			x <- summarise_glm(temp, j, i)
		} else {
			temp <- lm(fom, data = dat)
			x <- summarise_lm(temp, j, i)
			x$summary_dat <- x$summary_dat %>%
				dplyr::select(-adj_r2)
		}
		comp_res[[i]] <- x
		sum_res[[i]] <- x$summary_dat
	}
	temp <- do.call(rbind, comp_res)
	temp <- as.data.frame(temp) %>%
		rownames_to_column(var = "snp") %>%
		mutate(snp = gsub("_[ACTG]", "", snp))
	comp_res2[[j]] <- temp

	temp <- do.call(rbind, sum_res)
	temp <- as.data.frame(temp) %>%
		rownames_to_column(var = "snp") %>%
		mutate(effect_allele = gsub(".*_", "", snp)) %>%
		mutate(snp = gsub("_[ACTG]", "", snp))
	sum_res2[[j]] <- temp
}

res <- do.call(rbind, sum_res2)
rownames(res) <- NULL

write.table(res, file = paste0("data/", paste(traits, collapse = "_"), "_spouse_diff_gwas_res.txt"), col.names = T, row.names = F, qu = F, sep = "\t")

