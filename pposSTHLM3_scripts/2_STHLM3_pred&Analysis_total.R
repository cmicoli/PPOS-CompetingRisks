#==============================================================================
# FILENAME: 2_STHLM3_pred&Analysis_total.R
# PROJECT: 	STHLM3 PPOS
# PURPOSE: Predicting & analysis different model strategies (PHASE 2 + 3 of PPOS) to data (PCH and weibull models)
# AUTHOR: Chiara Micoli

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)
#==============================================================================

## PHASE 2/3 of PPOS: Predicting times after fitting models to STHLM3 data, 
## and analyzing the results of the prediction
## For each modelling strategy, the results of the prediction are saved in a list with 
## the time added and the results of the prediction

############################# model strategy A  - diff splits ##########################
rm(list = ls())
print("model strategy A - diff splits")
load("/nfs/home/chimic/STHLM3/fit_pch_sthlm3_diffSplits_stan.RData")
splits <- "diff"
tchange1 <- c(0, tsplit1)
tchange2 <- c(0, tsplit2)
savedir <-  "/nfs/home/chimic/STHLM3/PRED_rep_sthlm3_diffSplits_stan.RData" 
source("/nfs/home/chimic/STHLM3/2a_STHLM3_pred&analysis_pch.R", echo=TRUE)

############################# model strategy B - 2 splits ##########################
rm(list = ls())
print("model strategy B - 2 splits")
load("/nfs/home/chimic/STHLM3/fit_pch_sthlm3_2Splits_stan.RData")
splits <- 2
tchange1 <- c(0, tsplit1)
tchange2 <- c(0, tsplit2)
savedir <- "/nfs/home/chimic/STHLM3/PRED_rep_sthlm3_2Splits_stan.RData"
source("/nfs/home/chimic/STHLM3/2a_STHLM3_pred&analysis_pch.R", echo=TRUE)

############################# model strategy C - 4 splits ##########################
rm(list = ls())
print("model strategy C - 4 splits")
load("/nfs/home/chimic/STHLM3/fit_pch_sthlm3_4Splits_stan.RData")
splits <- 4 
tchange1 <- c(0, tsplit1)
tchange2 <- c(0, tsplit2)
savedir <- "/nfs/home/chimic/STHLM3/PRED_rep_sthlm3_4Splits_stan.RData"
source("/nfs/home/chimic/STHLM3/2a_STHLM3_pred&analysis_pch.R", echo=TRUE)

############################# model strategy D - weibull ##########################
rm(list = ls())
print("model strategy D - weibull")
load("/nfs/home/chimic/STHLM3/fit_weibull_sthlm3_1016.RData")
savedir <- "/nfs/home/chimic/STHLM3/PRED_rep_weibull_sthlm3_1016.RData"
source("/nfs/home/chimic/STHLM3/2b_STHLM3_pred&analysis_weibull.R", echo=TRUE) 



sessionInfo()