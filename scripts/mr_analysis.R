# ------------------------------------------------
# Within spouse MR analyses
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "data.table", "TwoSampleMR", "MRInstruments")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

# ------------------------------------------------
# read in exposure and outcome data
# ------------------------------------------------
traits <- c("chd", "height")
outcome <- "chd"

dat <- read_delim("data/chd_height_spouse_diff_gwas_res.txt", delim = "\t")
dat <- read_delim("data/differences_dat_chd_height.txt", delim = "\t")

exp_dat <- dat %>%
	dplyr::filter(outcome == "height") %>%
	dplyr::filter(p < 0.05)

exp_dat <- format_data(
	dat = exp_dat,
	type = "exposure",
	phenotype_col = "outcome",
	snp_col = "snp", 
	beta_col = "estimate", 
	se_col = "se", 
	pval_col = "p", 
	effect_allele_col = "effect_allele", 
	other_allele_col = "T"
)
head(exp_dat)

out_dat <- dat %>%
	dplyr::filter(snp %in% exp_dat$snp) %>%
	dplyr::filter(outcome == "chd")

out_dat <- format_data(
	dat = out_dat,
	type = "outcome",
	phenotype_col = "outcome",
	snp_col = "snp", 
	beta_col = "estimate", 
	se_col = "se", 
	pval_col = "p", 
	effect_allele_col = "effect_allele",
	other_allele_col = "T"
)

dat <- harmonise_data(exp_dat, out_dat, action = 1)
res <- mr(dat)

# ------------------------------------------------
# run mr 
# ------------------------------------------------
g0 = g[y==0]
tsps.glm = glm(y~predict(lm(x[y==0]~g0), newdata=list(g0=g)),
family=binomial)
beta_tsps = tsps.glm$coef[2]
se_tsps = summary(tsps.glm)$coef[2,2]

snps <- grep("rs[0-9]", colnames(dat), value = T)
snp_list <- paste(snps, collapse = "+")
no_out <- dat %>%
	dplyr::filter(chd_diff == 0)

tsps.glm <- glm(
	dat$chd_diff ~ predict(lm(as.formula(paste0("height_diff ~ ", snp_list)), data = no_out), newdata = list(no_out=dat)),
	family = binomial
	)

library(sem)
beta_tsls = tsls(y, cbind(x, rep(1,N)), cbind(g, rep(1,N)),w=rep(1,N))$coef[1]
se_tsls = sqrt(tsls(y, cbind(x, rep(1,N)), cbind(g, rep(1,N)),
w=rep(1,N))$V[1,1])
library(ivpack)
ivmodel = ivreg(y~x|g, x=TRUE)
summary(ivmodel)
beta_tsls = ivreg(y~x|g, x=TRUE)$coef[2]
se_tsls = summary(ivreg(y~x|g, x=TRUE))$coef[2,2]

mr_wald <- out_dat$estimate / exp_dat$estimate

pdf("data/output/test.pdf")
hist(mr_wald)
dev.off()

# ------------------------------------------------
# read in phenos from mr base
# ------------------------------------------------
# needs to be done not on bc because mr fucking base doesn't fucking work on here the fuck
load("data/giant_height_snps.RData")
g_height <- x %>%
	dplyr::filter(pval.exposure < 1e-8) %>%
	clump_data()
rm(x)

ao <- available_outcomes()
ao[grep(paste(c("coronary"), collapse = "|"), ao$trait, ignore.case = T),]
out <- extract_outcome_data(snps = g_height$SNP, outcomes = 7)
head(out)
dat <- harmonise_data(
    exposure_dat = g_height, 
    outcome_dat = out
)
head(dat)
res <- mr(dat)

# 1. Read in both exposure and outcome gwas dat
# 2. Run the MR analyses - sensitivity analyses? 
# 3. Read in phenos from MR base
# 4. Estimates from step 2 with estimates from step 3
# 5. Party 