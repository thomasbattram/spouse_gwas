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
# dat <- read_delim("data/differences_dat_chd_height.txt", delim = "\t")

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
	dplyr::filter(snp %in% exp_dat$SNP) %>%
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

# ------------------------------------------------
# run mr 
# ------------------------------------------------

dat <- harmonise_data(exp_dat, out_dat, action = 1)
res <- mr(dat)
scat <- mr_scatter_plot(res, dat)
ggsave("data/output/mr_chd_height_scatter.pdf", plot = scat[[1]])

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