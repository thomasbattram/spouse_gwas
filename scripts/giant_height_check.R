# ------------------------------------------------
# running the analysis for spouse gwas stuff
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "ggExtra", "TwoSampleMR")
lapply(pkgs, require, character.only = T)


########
# Do bit below out of bluecrystal!!
# ao <- available_outcomes()
# head(ao)

# ao[grep("height", ao$trait, ignore.case = T), ]
# 25282103

# x <- extract_instruments(outcomes = 89, clump = F)
# save(x, file = "~/Desktop/spouse_gwas/giant_height_snps.RData")
########
devtools::load_all("~/repos/usefunc/")

# 1. get data from gwas catalog
load("data/giant_height_snps.RData")
g_height <- x
colnames(g_height)
g_height <- g_height %>%
	mutate(snp = SNP) %>%
	mutate(estimate = beta.exposure) %>%
	mutate(se = se.exposure) %>%
	mutate(p = pval.exposure) %>%
	mutate(a1 = effect_allele.exposure) %>%
	mutate(cohort = "giant")

# 2. compare ukb conventional gwas res with giant res
ukb_height <- read_table("data/output/Height-Height.conventional.Height.assoc.linear")
ukb_height <- ukb_height %>%
	mutate(snp = SNP) %>%
	mutate(estimate = BETA) %>%
	mutate(se = SE) %>%
	mutate(p = P) %>% 
	mutate(cohort = "ukb") %>%
	mutate(a1 = A1) %>%
	separate(`L95      U95`, into = c("CI_low", "CI_up"), sep = "\\s+")

temp_ukb <- ukb_height %>%
	dplyr::select(snp, estimate, se, p, a1, cohort) %>%
	dplyr::filter(snp %in% g_height$snp)

temp_g <- g_height %>%
	dplyr::select(snp, estimate, se, p, a1, cohort) %>%
	dplyr::filter(snp %in% temp_ukb$snp)

# allele check
m_dat <- temp_g %>%
	left_join(temp_ukb, by = c("snp" = "snp"))

al_diff <- m_dat[m_dat$a1.x != m_dat$a1.y, "snp"]
dplyr::filter(temp_g, snp %in% al_diff)
for (i in m_dat$snp) {
	if (i %in% al_diff) {
		print(i)
		temp_g[temp_g$snp == i, "estimate"] <- temp_g[temp_g$snp == i, "estimate"] * -1
	}
}
dplyr::filter(temp_g, snp %in% al_diff)

fin_dat <- rbind(temp_g, temp_ukb) %>%
	arrange(snp)

dplyr::filter(fin_dat, snp %in% al_diff)

ukb_est <- fin_dat[fin_dat$cohort == "ukb", "estimate"]
giant_est <- fin_dat[fin_dat$cohort == "giant", "estimate"]

est_diff <- ukb_est - giant_est

hdat <- data.frame(snp = unique(fin_dat$snp), ukb_est, giant_est, est_diff)

rm(list = c("ukb_est", "giant_est", "est_diff"))

fit <- lm(giant_est ~ ukb_est, data = hdat)
summary(fit)
splot <- gglot(hdat, aes(x = ukb_est, y = giant_est)) +
	geom_point() +
	# geom_abline(intercept = 0, slope = 1, colour = "red") +
	geom_smooth(method = "lm")

ggsave("data/output/giant_ukb_comparison.pdf", ggMarginal(splot, type = "histogram", xparams = list(bins = 50), yparams = list(bins = 50)))


# - If difference in effect estimates is ~0 move onto 3
# 3. compare within spouse differences with giant height snps


