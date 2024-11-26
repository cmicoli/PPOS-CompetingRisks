#==============================================================================
# FILENAME: 1_main_noPH.R
# PROJECT: 	ISPY-COVID PPOS
# PURPOSE: Main analysis for PPOS estimation, stratifying modeling over treatment arm
# AUTHOR: Chiara Micoli, Andrea Discacciati

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)
#==============================================================================
## Main analysis for PPOS estimation, stratifying modeling over treatment arm. 
## The analysis is done for the ISPY-COVID data, 
# and the main goal is to estimate the PPOS based on the observed data. 
## PHASE 1/2/3 are all done in this script. Final result is the ppos. 

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
library(survival)
library(parallel)
library(rms)
library(MASS)
library(brms)
library(dplyr)
library(purrr)
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")

### parameters to change
cores <- 40
nsim <- 2500
savedir <- "/nfs/home/chimic/ISPYCOVID/fit_main_A2500_20241106.RData"

## loading functions for analysis
source("/nfs/home/chimic/ISPYCOVID/00Functions_weibull_sim.R")
## loading data
cyclosporine.timemachine <- readRDS("/nfs/home/chimic/ISPYCOVID/20240201_cyclosporine.timemachine.RDS")


################################################################################
# splitting data in d_obs and d_cens
dat <- cyclosporine.timemachine
dat$EVENT <- factor(dat$rec_EV, 
                    levels = c("Censored", "Recovery", "Death w/o recovery"), 
                    labels = c(0,1,2))
dat$TREATMENT <- factor(dat$ARM,
                        levels = c("Control", "Cyclosporine"), 
                        labels = c("CNTR", "NEW_TR"))
dat$event.all <- ifelse(dat$rec_EV %in% c("Recovery", "Death w/o recovery"), 1, 0)
dat$event1 <- dat$event.all * (dat$rec_EV == "Recovery")
dat$event2 <- dat$event.all * (dat$rec_EV == "Death w/o recovery")
dat$censored1 <- 1 - dat$event1
dat$censored2 <- 1 - dat$event2
with(dat, prop.table(table(rec_EV)))

dat_new.tr <- dat %>% filter(TREATMENT == "NEW_TR")
dat_control <- dat %>% filter(TREATMENT == "CNTR")

# splitting between completely observed and censored data
dat1 <- dat %>% filter(event.all == 1 | rec_FU  == 60) %>% 
  select(USUBJID, TREATMENT, SCRNCOVID_bin, rec_FU, EVENT)  ## 101 patients
dat0_cens <- dat %>% filter(event.all == 0 & rec_FU <60) ## 32 patients

# predicting baseline for new data (patients that will be recruted - D_new)
## new treatment has to be present 67 times, X is the number of controls needed to enroll 67 new_tr
D_new <- lapply (1:nsim, function(x) {
  n.control <- rnbinom(1, 67, 0.45)
  d <- data.frame(
    USUBJID = paste0("NEW", seq_len(n.control+ 67)),
    TREATMENT = factor(rep(c("CNTR", "NEW_TR"), times = c(n.control, 67))),
    SCRNCOVID_bin = factor(sample(c("5", "6/7"), n.control+67, replace = TRUE, prob = c(0.9, 0.1)))
  )
})

# dat0 containes data that need prediction, either update of follow up for dat0_cens, or new pred for D_new
dat0 <- lapply (D_new, function(x) {
  dat0 <- dat0_cens %>% 
    select(USUBJID, TREATMENT, SCRNCOVID_bin, rec_FU, EVENT) %>%
    rbind (x %>% 
             mutate(rec_FU = 0, EVENT = 0))
})

# new treatment datasets
dat_new.tr1 <- dat_new.tr %>% filter(event.all == 1 | rec_FU  == 60) %>% 
  select(USUBJID, TREATMENT, SCRNCOVID_bin, rec_FU, EVENT)  ## 42 patients
dat_new.tr0 <- lapply(dat0, \(x) x %>% filter(TREATMENT == "NEW_TR") )
                      #function(x) {  dat_new.tr0 <-  x %>% filter(TREATMENT == "NEW_TR")})
## 83 patients to predict: 16 already enrolled, 67 new

# control datasets
dat_control1 <- dat_control %>% filter(event.all == 1 | rec_FU  == 60) %>% 
  select(USUBJID, TREATMENT, SCRNCOVID_bin, rec_FU, EVENT)  ## 59 patients
dat_control0 <- lapply(dat0, \(x) x %>% filter(TREATMENT == "CNTR") )
## XX patients to predict: 16 already enrolled, new patients in the trial (new accruals - from negative binomial)


#####################################################################################
# PHASE 1: modelling
## models specification
bbff1 <- bf(rec_FU | cens(censored1) ~   0 + Intercept  + SCRNCOVID_bin,
            family = weibullPH)
bbff2 <- bf(rec_FU | cens(censored2) ~   0 + Intercept  + SCRNCOVID_bin,
            family = weibullPH)

# model fit for new treatment 
fit_new.tr <- fit.noPH(dataset = dat_new.tr, 
                       brms_function1 = bbff1,
                       brms_function2 = bbff2   )

# model fit for control arm 
fit_control <- fit.noPH(dataset = dat_control, 
                        brms_function1 = bbff1,
                        brms_function2 = bbff2   )

#####################################################################################
# PHASE 2: prediction
## new predicted dataset 
predicted_new.tr <- new_predicted_dataset_noPH(data_to_predict = dat_new.tr0, 
                                               data_observed = dat_new.tr1, 
                                               n.sim = nsim, 
                                               model1 = fit_new.tr$Fit1, 
                                               model2 = fit_new.tr$Fit2, 
                                               cores)

predicted_control <- new_predicted_dataset_noPH(data_to_predict = dat_control0, 
                                                data_observed = dat_control1, 
                                                n.sim = nsim, 
                                                model1 = fit_control$Fit1, 
                                                model2 = fit_control$Fit2, 
                                                cores)
# merging predicted datasets
predicted <- lapply (1:nsim , \(x) rbind(predicted_control[[x]], 
                                         predicted_new.tr [[x]]))

#####################################################################################
# PHASE 3: analysis
## analysis of new predicted dataset 
fit.null_analysis <- compile.analysis.model(predicted[[1]])
predicted_an <- mclapply(predicted, 
                         \(x)  main.analysis.update(dataset = x, fit.null = fit.null_analysis),
                         mc.cores = cores)

analysis.HR <- as.data.frame(do.call(rbind, predicted_an %>% map(2)))

## PPOS
table(analysis.HR$success)
table(analysis.HR$pr>0.975) %>% prop.table() 

#####################################################################################
# saving results
save.image(savedir)
