#==============================================================================
# FILENAME: 1a_STHLM3_modelFit_PCH_AR1.R
# PROJECT: 	STHLM3 PPOS
# PURPOSE: fitting pch models to sthlm3 data
# AUTHOR: Chiara Micoli, Andrea Discacciati

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)

#==============================================================================
## pch fit using stan with different or equal time splits for pca hazard and other cause hazard
# used for model strategy A,B,C of the paper

# libraries 
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

library(survival)
library(parallel)
library(rms)
library(cmdstanr)
library(MASS)
library(brms)
library(dplyr)
library(flexsurv)
library(tidybayes)
library(tidyverse)
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
PCH_AR1_STHLM3 <- cmdstan_model("/nfs/home/chimic/STHLM3/PCH_AR1_STHLM3.stan") 


### analysis done with different models based on invitation group ########
dat_tot <- person_main_surv[, c("LOPNR", "group_rand", "age_at_invitation", "date_invitation",
                                "time18", "event18", "event18.all" , "maxFU18"  , "Emigration_date"       )]
dat_tot$age_at_invitation_c <- dat_tot$age_at_invitation - mean(dat_tot$age_at_invitation)
summary(as.numeric(dat_tot$time18))
dat_inv <- dat_tot %>% filter(group_rand == "invited")
dat_not.inv <- dat_tot %>% filter(group_rand == "notinvited")


############### splitting data based on tsplit2  -> for prostate cancer mortality ######################
simdata.split.invPCA <- survSplit(Surv(time18, event18.all) ~ ., dat_inv, cut = tsplit2)
simdata.split.invPCA$event1 <- simdata.split.invPCA$event18.all * (simdata.split.invPCA$event18 == "1")
simdata.split.invPCA$event2 <- simdata.split.invPCA$event18.all * (simdata.split.invPCA$event18 == "2")
simdata.split.invPCA$censored1 <- 1 - simdata.split.invPCA$event1
simdata.split.invPCA$censored2 <- 1 - simdata.split.invPCA$event2
simdata.split.invPCA$tstart.c <- as.factor(simdata.split.invPCA$tstart)
with(dat_inv, prop.table(table(event18)))


simdata.split.not.invPCA <- survSplit(Surv(time18, event18.all) ~ ., dat_not.inv, cut = tsplit2)
simdata.split.not.invPCA$event1 <- simdata.split.not.invPCA$event18.all * (simdata.split.not.invPCA$event18 == "1")
simdata.split.not.invPCA$event2 <- simdata.split.not.invPCA$event18.all * (simdata.split.not.invPCA$event18 == "2")
simdata.split.not.invPCA$censored1 <- 1 - simdata.split.not.invPCA$event1
simdata.split.not.invPCA$censored2 <- 1 - simdata.split.not.invPCA$event2
simdata.split.not.invPCA$tstart.c <- as.factor(simdata.split.not.invPCA$tstart)
with(dat_not.inv, prop.table(table(event18)))

############### splitting data based on tsplit1  -> for other cause mortality ######################
simdata.split.invOC <- survSplit(Surv(time18, event18.all) ~ ., dat_inv, cut = tsplit1)
simdata.split.invOC$event1 <- simdata.split.invOC$event18.all * (simdata.split.invOC$event18 == "1")
simdata.split.invOC$event2 <- simdata.split.invOC$event18.all * (simdata.split.invOC$event18 == "2")
simdata.split.invOC$censored1 <- 1 - simdata.split.invOC$event1
simdata.split.invOC$censored2 <- 1 - simdata.split.invOC$event2
simdata.split.invOC$tstart.c <- as.factor(simdata.split.invOC$tstart)
#contrasts(simdata.split.invOC$tstart.c) <- codingMatrices::contr.diff(12)  # MASS::contr.sdif(number or ) 
with(dat_inv, prop.table(table(event18)))

simdata.split.not.invOC <- survSplit(Surv(time18, event18.all) ~ ., dat_not.inv, cut = tsplit1)
simdata.split.not.invOC$event1 <- simdata.split.not.invOC$event18.all * (simdata.split.not.invOC$event18 == "1")
simdata.split.not.invOC$event2 <- simdata.split.not.invOC$event18.all * (simdata.split.not.invOC$event18 == "2")
simdata.split.not.invOC$censored1 <- 1 - simdata.split.not.invOC$event1
simdata.split.not.invOC$censored2 <- 1 - simdata.split.not.invOC$event2
simdata.split.not.invOC$tstart.c <- as.factor(simdata.split.not.invOC$tstart)
#contrasts(simdata.split.not.invOC$tstart.c) <- codingMatrices::contr.diff(12) 
with(dat_not.inv, prop.table(table(event18)))


## how are the events splited
table(person_main_surv$event18, person_main_surv$group_rand)
table(simdata.split.invPCA$event1, simdata.split.invPCA$tstart)
table(simdata.split.invPCA$event2, simdata.split.invPCA$tstart)

table(simdata.split.not.invPCA$event1, simdata.split.not.invPCA$tstart)
table(simdata.split.not.invPCA$event2, simdata.split.not.invPCA$tstart)

table(simdata.split.not.invOC$event1, simdata.split.not.invOC$tstart)
table(simdata.split.invOC$event1, simdata.split.invOC$tstart)

# observed incidence rate
## event 1
ir.ni <- tapply(simdata.split.not.invOC$event1, simdata.split.not.invOC$tstart.c, sum) /
  tapply(simdata.split.not.invOC$time18-simdata.split.not.invOC$tstart, simdata.split.not.invOC$tstart.c, sum)
log(sapply(ir.ni, `[`, 1))
ir.i <- tapply(simdata.split.invOC$event1, simdata.split.invOC$tstart.c, sum) /
  tapply(simdata.split.invOC$time18-simdata.split.invOC$tstart, simdata.split.invOC$tstart.c, sum)
log(sapply(ir.i, `[`, 1))

## event 2
ir.ni <- tapply(simdata.split.not.invPCA$event2, simdata.split.not.invPCA$tstart.c, sum) /
  tapply(simdata.split.not.invPCA$time18-simdata.split.not.invPCA$tstart, simdata.split.not.invPCA$tstart.c, sum)
log(sapply(ir.ni, `[`, 1))
ir.i <- tapply(simdata.split.invPCA$event2, simdata.split.invPCA$tstart.c, sum) /
  tapply(simdata.split.invPCA$time18-simdata.split.invPCA$tstart, simdata.split.invPCA$tstart.c, sum)
log(sapply(ir.i, `[`, 1))

##--> use the observed log(ir) to center the prior


# Design matrix for invited men (for the two cause of failure)
X.invPCA <- model.matrix(~ -1 + tstart.c + age_at_invitation_c, data = simdata.split.invPCA)
stan.simdata.split.invPCA <- list( # Prepare data for Stan
  N = nrow(simdata.split.invPCA),
  Y = simdata.split.invPCA$time18,
  cens = simdata.split.invPCA$censored2,
  lb = simdata.split.invPCA$tstart,
  K = ncol(X.invPCA),
  X = X.invPCA,
  prior_only = 0
)
X.invOC <- model.matrix(~ -1 + tstart.c + age_at_invitation_c, data = simdata.split.invOC)
stan.simdata.split.invOC <- list(  # Prepare data for Stan
  N = nrow(simdata.split.invOC),
  Y = simdata.split.invOC$time18,
  cens = simdata.split.invOC$censored1,
  lb = simdata.split.invOC$tstart,
  K = ncol(X.invOC),
  X = X.invOC,
  prior_only = 0
)

# Design matrix for non-invited men (for the two cause of failure)
X.not.invPCA <- model.matrix(~ -1 + tstart.c + age_at_invitation_c, data = simdata.split.not.invPCA)
stan.simdata.split.not.invPCA <- list( # Prepare data for Stan
  N = nrow(simdata.split.not.invPCA),
  Y = simdata.split.not.invPCA$time18,
  cens = simdata.split.not.invPCA$censored2,
  lb = simdata.split.not.invPCA$tstart,
  K = ncol(X.not.invPCA),
  X = X.not.invPCA,
  prior_only = 0
)
X.not.invOC <- model.matrix(~ -1 + tstart.c + age_at_invitation_c, data = simdata.split.not.invOC)
stan.simdata.split.not.invOC <- list(  # Prepare data for Stan
  N = nrow(simdata.split.not.invOC),
  Y = simdata.split.not.invOC$time18,
  cens = simdata.split.not.invOC$censored1,
  lb = simdata.split.not.invOC$tstart,
  K = ncol(X.not.invOC),
  X = X.not.invOC,
  prior_only = 0
)


#function to save models - from stan
cmdstanr.resolve <- function(fit){
  temp_rds_file <- tempfile(fileext = ".RDS")
  fit$save_object(file = temp_rds_file)
  fit <- readRDS(temp_rds_file)  
  return(fit)
}

## fit OTHER CAUSE mortality - NOT INVITED
fit1_not.inv <- PCH_AR1_STHLM3$sample(data = stan.simdata.split.not.invOC, 
                                      refresh = 0,
                                      chains = CHAINS,
                                      seed = SEED,
                                      iter_warmup = ITER.WARMUP,
                                      iter_sampling = ITER.SAMPLING,
                                      show_exceptions = FALSE)
fit1_not.inv <- cmdstanr.resolve(fit1_not.inv)
fit1_not.inv_fsr <- flexsurvreg(Surv(tstart, time18, event1) ~ tstart.c + age_at_invitation_c,
                                data = simdata.split.not.invOC,
                                dist = "exp")



## fit PCA mortality - NOT INVITED

fit2_not.inv <- PCH_AR1_STHLM3$sample(data = stan.simdata.split.not.invPCA, 
                                      refresh = 0,
                                      chains = CHAINS,
                                      seed = SEED,
                                      iter_warmup = ITER.WARMUP,
                                      iter_sampling = ITER.SAMPLING,
                                      show_exceptions = FALSE)
fit2_not.inv <- cmdstanr.resolve(fit2_not.inv)
fit2_not.inv_fsr <- flexsurvreg(Surv(tstart, time18, event2) ~ tstart.c + age_at_invitation_c,
                                data = simdata.split.not.invPCA,
                                dist = "exp")



## fit OTHER CAUSE mortality - INVITED

fit1_inv <- PCH_AR1_STHLM3$sample(data = stan.simdata.split.invOC, 
                                  refresh = 0,
                                  chains = CHAINS,
                                  seed = SEED,
                                  iter_warmup = ITER.WARMUP,
                                  iter_sampling = ITER.SAMPLING,
                                  show_exceptions = FALSE)
fit1_inv <- cmdstanr.resolve(fit1_inv)
fit1_inv_fsr <- flexsurvreg(Surv(tstart, time18, event1) ~ tstart.c + age_at_invitation_c,
                            data = simdata.split.invOC,
                            dist = "exp")

#------

## fit PCA mortality - INVITED

fit2_inv <- PCH_AR1_STHLM3$sample(data = stan.simdata.split.invPCA, 
                                  refresh = 0,
                                  chains = CHAINS,
                                  seed = SEED,
                                  iter_warmup = ITER.WARMUP,
                                  iter_sampling = ITER.SAMPLING,
                                  show_exceptions = FALSE)
fit2_inv <- cmdstanr.resolve(fit2_inv)
fit2_inv_fsr <- flexsurvreg(Surv(tstart, time18, event2) ~ tstart.c + age_at_invitation_c,
                            data = simdata.split.invPCA,
                            dist = "exp")


## saving results
save.image(savedir)

