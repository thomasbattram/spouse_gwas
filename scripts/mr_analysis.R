# ------------------------------------------------
# Within spouse MR analyses
# ------------------------------------------------

rm(list = ls())

setwd("~/spouse_gwas")

pkgs <- c("tidyverse", "conflicted", "data.table")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

# ------------------------------------------------
# read in exposure and outcome data
# ------------------------------------------------
exposure <- "height"
outcome <- "chd"

exp_dat <- read_delim(paste0(exposure, "_spouse_diff_gwas_res.txt"), delim = "\t")
out_dat <- read_delim(paste0(outcome, "_spouse_diff_gwas_res.txt"), delim = "\t")

# extract replicated snps
exp_dat <- exp_dat %>%
	dplyr::filter(p < 0.05)

out_dat <- out_dat %>%
	dplyr::filter(snp %in% exp_dat$snp)

mr_wald <- out_dat$estimate / exp_dat$estimate



# 1. Read in both exposure and outcome gwas dat
# 2. Run the MR analyses - sensitivity analyses? 
# 3. Read in phenos from MR base
# 4. Estimates from step 2 with estimates from step 3
# 5. Party 