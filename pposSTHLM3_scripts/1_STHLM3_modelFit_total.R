#==============================================================================
# FILENAME: 1_STHLM3_modelFit_total.R
# PROJECT: 	STHLM3 PPOS
# PURPOSE: Fitting different model strategies (PHASE 1 of PPOS) to data (PCH and weibull models)
# AUTHOR: Chiara Micoli

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)

#==============================================================================
## doing different modeling strategy for PHASE 1 of PPOS and saving them


# Compile Stan program (change path as needed)
CHAINS = 5
SEED = 12345
ITER.WARMUP = 500
ITER.SAMPLING = 1000
CORES = 5 

load(file = "/nfs/home/chimic/STHLM3/person_main_surv_KOSMOS_241022.RData")
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
library(brms)

################## parameters that change in each modelling strategy #################

############################# model strategy A ##########################
print("fit model strategy A - diff splits")
tsplit1 <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5) ## splits for other cause mortality
tsplit2 <- c(2,3,4,5) ## splits for pca mortality
savedir <- "/nfs/home/chimic/STHLM3/fit_pch_sthlm3_diffSplits_stan.RData"
source("/nfs/home/chimic/STHLM3/1a_STHLM3_modelFit_PCH_AR1.R", echo=TRUE)

############################# model strategy B ##########################
print("fit model strategy B - 2 splits")
tsplit1 <- c(2,4) ## splits for other cause mortality
tsplit2 <- c(2,4) ## splits for pca mortality
savedir <- "/nfs/home/chimic/STHLM3/fit_pch_sthlm3_2Splits_stan.RData"
source("/nfs/home/chimic/STHLM3/1a_STHLM3_modelFit_PCH_AR1.R", echo=TRUE)

############################# model strategy C ##########################
print("fit model strategy C - 4 splits")
tsplit1 <- c(2,3,4,5) ## splits for other cause mortality
tsplit2 <- c(2,3,4,5) ## splits for pca mortality
savedir <- "/nfs/home/chimic/STHLM3/fit_pch_sthlm3_4Splits_stan.RData"
source("/nfs/home/chimic/STHLM3/1a_STHLM3_modelFit_PCH_AR1.R", echo=TRUE)

############################# model strategy D -- Weibull ##########################
print("fit model strategy D - weibull")
DISTRIBUTION <- weibullPH
DISTRIBUTION.name <- "weibullPH"
priors1 <-  priors2 <- set_prior("normal(-10, 20)", class =  "b", coef = "Intercept") +                  ## scale --> log(0.004713   )
  set_prior("normal(0, sqrt(0.5))", class = "b", coef = "age_at_invitation_c") +
  set_prior("exponential(1)", class = "gamma") 

savedir <- "/nfs/home/chimic/STHLM3/fit_weibull_sthlm3_1016.RData"
source("/nfs/home/chimic/STHLM3/1b_STHL3_modelFit_brms.R", echo=TRUE)




sessionInfo()