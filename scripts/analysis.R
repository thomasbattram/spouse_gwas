# ------------------------------------------------
# running the analysis for spouse gwas stuff
# ------------------------------------------------
rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "ggExtra")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")
trait <- "height"
dat <- read_delim(paste0("data/differences_dat_", trait, ".txt"), delim = "\t")

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
	fom <- as.formula(paste0("trait_diff ~ ", i, " + age_diff + Sex"))
	temp <- lm(fom, data = dat)
	x <- summarise_lm(temp, trait, i)
	comp_res[[i]] <- x
	sum_res[[i]] <- x$summary_dat	
}

res <- do.call(rbind, sum_res)
res <- as.data.frame(res) %>%
	rownames_to_column(var = "snp") %>%
	mutate(snp = gsub("_[ACTG]", "", snp))

write.table(res, file = paste0("data/", trait, "_spouse_diff_gwas_res.txt"), col.names = T, row.names = F, qu = F, sep = "\t")

################################################
##
##   MAKE BIT BELOW INTO DIFFERENT SCRIPT!!! ### also make bit above different - hard to get comparison and mr scripts the same
##
################################################



res_sig <- res %>%
	dplyr::filter(p < 0.05) %>%
	dplyr::select(-outcome, -adj_r2) %>%
	mutate(mod = "spouse_only")

# read in actual gwas for comparison
height_gwas <- read_table("data/output/Height-Height.conventional.Height.assoc.linear")
g_dat <- height_gwas %>%
	mutate(snp = SNP) %>%
	mutate(estimate = BETA) %>%
	mutate(se = SE) %>%
	mutate(p = P) %>% 
	mutate(mod = "conventional") %>%
	separate(`L95      U95`, into = c("CI_low", "CI_up"), sep = "\\s+") %>%
	dplyr::select(snp, estimate, se, p, CI_low, CI_up, mod) %>%
	dplyr::filter(snp %in% res_sig$snp)

# Allele check
sp_a1 <- gsub(".*_", "", snps)
sp_snp <- gsub("_.*", "", snps)
temp <- data.frame(snp = sp_snp, sp_a1)
al_check <- temp %>%
	left_join(height_gwas, by = c("snp" = "SNP")) %>%
	dplyr::select(snp, sp_a1, A1) 

al_diff <- al_check[al_check$sp_a1 != al_check$A1, "snp"]

dplyr::filter(g_dat, snp %in% al_diff)
for (i in g_dat$snp) {
	if (i %in% al_diff) {
		print(i)
		g_dat[g_dat$snp == i, "estimate"] <- g_dat[g_dat$snp == i, "estimate"] * -1
	}
}
dplyr::filter(g_dat, snp %in% al_diff)

# join data together - probs do it so they're side by side then remove everything but beta, se, p, l95, u95 for both 

fin_dat <- rbind(res_sig, g_dat)
head(fin_dat)

fin_dat <- fin_dat %>%
	mutate(facet_var = facet_var_gen(fin_dat, 6, "mod")) %>%
	mutate(mod = as.factor(fin_dat$mod)) %>%
	mutate(`2.5 %` = as.numeric(CI_low)) %>%
	mutate(`97.5 %` = as.numeric(CI_up)) %>%
	mutate(Estimate = estimate) %>%
	dplyr::select(-estimate, -CI_up, -CI_low) %>%
	arrange(snp)

formals(forest_plot)

fplot <- forest_plot(fin_dat, col_num = 6, group = "mod", y_axis = "snp", null_at = 0)

ggsave("data/output/height_comp.pdf", plot = fplot, width = 15, height = 10, units = "in")

spo_est <- fin_dat[fin_dat$mod == "spouse_only", "Estimate"]
con_est <- fin_dat[fin_dat$mod == "conventional", "Estimate"]

est_diff <- con_est - spo_est

hdat <- data.frame(snp = unique(fin_dat$snp), spo_est, con_est, est_diff)


rm(list = c("spo_est", "con_est", "est_diff"))
hplot <- ggplot(hdat) +
	geom_histogram(aes(x = est_diff), alpha = 0.2) +
	geom_histogram(aes(x = spo_est), fill = "red", alpha = 0.2) +
	geom_histogram(aes(x = con_est), fill = "blue", alpha = 0.2)

ggsave("data/output/height_diff.pdf", plot = hplot, width = 15, height = 10, units = "in")

fit <- lm(spo_est ~ con_est, data = hdat)
summary(fit)
splot <- ggplot(hdat, aes(x = con_est, y = spo_est)) +
	geom_point() +
	geom_abline(intercept = 0, slope = 1, colour = "red") +
	geom_smooth(method = "lm")

ggsave("data/output/height_diff_scat.pdf", ggMarginal(splot, type = "histogram", xparams = list(bins = 50), yparams = list(bins = 50)))


# fin_dat <- res_sig %>%
# 	left_join(height_gwas, by = c("snp" = "SNP"))


