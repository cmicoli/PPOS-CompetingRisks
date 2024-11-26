#==============================================================================
# FILENAME: 1_STHLM3_modelFit_brms.R
# PROJECT: 	STHLM3 PPOS
# PURPOSE: fitting brms models to sthlm3 data
# AUTHOR: Chiara Micoli, Andrea Discacciati

# R VERSION: R version 4.3.2 (2023-10-31)
    #Platform: x86_64-pc-linux-gnu (64-bit)
    #Running under: Red Hat Enterprise Linux 9.4 (Plow)

#==============================================================================
## fit using brms for bayesian models (weibull, exponential), for pca hazard and other cause hazard
# used for model strategy D of the paper

# libraries 
RNGkind("L'Ecuyer-CMRG")
set.seed(124)

library(survival)
library(parallel)
library(rms)
library(MASS)
library(brms)
library(dplyr)
library(writexl)
library(rstanarm)
library(survival)
library(tidybayes)
library(tidyverse)
library(cmdstanr)
library(flexsurv)
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")



### analysis done with different models based on invitation group ########
dat_tot <- person_main_surv[, c("LOPNR", "group_rand", "age_at_invitation", "date_invitation",
                                "time18", "event18", "event18.all" , "maxFU18"  , "Emigration_date"       )]
dat_tot$age_at_invitation_c <- dat_tot$age_at_invitation - mean(dat_tot$age_at_invitation)
summary(as.numeric(dat_tot$time18))
# dat_tot <- dat_tot[sample(1000), ]
dat_tot$event18.all <- ifelse(dat_tot$event18 %in% c("1", "2"), 1, 0)
dat_tot$event1 <- dat_tot$event18.all * (dat_tot$event18 == "1")
dat_tot$event2 <- dat_tot$event18.all * (dat_tot$event18 == "2")
dat_tot$censored1 <- 1 - dat_tot$event1
dat_tot$censored2 <- 1 - dat_tot$event2
dat_tot$time18 <- as.numeric(dat_tot$time18)
dat_inv <- dat_tot %>% filter(group_rand == "invited")
dat_not.inv <- dat_tot %>% filter(group_rand == "notinvited")




## flexsurv reg fit with weibull
fit1_inv_fsr <- flexsurv::flexsurvreg(Surv( time18, event18 == "1") ~  age_at_invitation_c,
                                          data = dat_inv,
                                          dist = DISTRIBUTION.name)
fit2_inv_fsr <- flexsurv::flexsurvreg(Surv( time18, event18 == "2") ~  age_at_invitation_c,
                                          data = dat_inv,
                                          dist = DISTRIBUTION.name)

fit1_not.inv_fsr <- flexsurv::flexsurvreg(Surv( time18, event18 == "1") ~  age_at_invitation_c,
                                              data = dat_not.inv,
                                              dist = DISTRIBUTION.name)

fit2_not.inv_fsr <- flexsurv::flexsurvreg(Surv( time18, event18 == "2") ~  age_at_invitation_c,
                                              data = dat_not.inv,
                                              dist = DISTRIBUTION.name)



## brms fit 
bbff1 <- bf(time18 | cens(censored1) ~   0 + Intercept + age_at_invitation_c,
            family = DISTRIBUTION)
bbff2 <- bf(time18 | cens(censored2) ~   0 + Intercept + age_at_invitation_c,
            family = DISTRIBUTION)

# priors1 <- set_prior("normal(-10, 20)", class =  "b", coef = "Intercept") +                  ## scale --> log(0.004713   )
#   set_prior("normal(0, sqrt(0.5))", class = "b", coef = "age_at_invitation_c") +
#   set_prior("exponential(1)", class = "gamma") 
# 
# priors2 <- set_prior("normal(-10, 20)", class =  "b", coef = "Intercept") +                 ## scale --> log(0.004713   )
#   set_prior("normal(0, sqrt(0.5))", class = "b", coef = "age_at_invitation_c") +
#   set_prior("exponential(1)", class = "gamma") 



fit1_inv <- brm(bbff1, data = dat_inv,     
                    chains =  CHAINS,
                    iter = ITER.SAMPLING + ITER.WARMUP,  ##1500
                    warmup = ITER.WARMUP, ## 500
                    seed = SEED,
                    cores = CORES,
                    stanvars = stanvars_weibullPH,
                    prior = priors1 ,
                    backend = "cmdstanr")

fit2_inv <- brm(bbff2, data = dat_inv,     
                    chains =  CHAINS,
                    iter = ITER.SAMPLING + ITER.WARMUP,  ##1500
                    warmup = ITER.WARMUP, ## 500
                    seed = SEED,
                    cores = CORES,
                    stanvars = stanvars_weibullPH,
                    prior = priors2 ,
                    backend = "cmdstanr")



fit1_not.inv <- brm(bbff1, data = dat_not.inv,     
                        chains =  CHAINS,
                        iter = ITER.SAMPLING + ITER.WARMUP,  ##1500
                        warmup = ITER.WARMUP, ## 500
                        seed = SEED,
                        cores = CORES,
                        stanvars = stanvars_weibullPH,
                        prior = priors1,
                        backend = "cmdstanr")

fit2_not.inv <- brm(bbff2, data = dat_not.inv,     
                        chains =  CHAINS,
                        iter = ITER.SAMPLING + ITER.WARMUP,  ##1500
                        warmup = ITER.WARMUP, ## 500
                        seed = SEED,
                        cores = CORES,
                        stanvars = stanvars_weibullPH,
                        prior = priors2,
                        backend = "cmdstanr")


save.image(savedir)