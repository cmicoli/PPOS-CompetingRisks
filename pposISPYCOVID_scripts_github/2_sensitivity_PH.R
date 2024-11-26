#==============================================================================
# FILENAME: 2_sensitivity_PH.R
# PROJECT: 	ISPY-COVID PPOS
# PURPOSE: Sensitivity PPOS analysis for the prior on the treatment effect in the weibullPH models
# AUTHOR: Chiara Micoli, Andrea Discacciati

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)
#==============================================================================
## Running everything with prior sd for treatment effect = 0.5 -> needs to be changed to 0.2, 0.1

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
sd.fit <-  0.5  ## 0.5, 0.1, 0.2
savedir <- "/nfs/home/chimic/ISPYCOVID/fit_predictions_priorPH_sd05_A2500_20241106.RData"
#"/nfs/home/chimic/ISPYCOVID/fit_predictions_priorPH_sd01_A2500_20241106.RData"
# "/nfs/home/chimic/ISPYCOVID/fit_predictions_priorPH_sd05_A2500_20241106.RData"
# "/nfs/home/chimic/ISPYCOVID/fit_predictions_priorPH_sd02_A2500_20241106.RData"

## loading functions for analysis
source("/nfs/home/chimic/ISPYCOVID/00Functions_weibull_sim.R")
## load data
cyclosporine.timemachine <- readRDS("/nfs/home/chimic/ISPYCOVID/20240201_cyclosporine.timemachine.RDS")

#####################################################################################
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


#####################################################################################
## models specification
bbff1 <- bf(rec_FU | cens(censored1) ~   0 + Intercept + TREATMENT + SCRNCOVID_bin,
            family = weibullPH)
bbff2 <- bf(rec_FU | cens(censored2) ~   0 + Intercept + TREATMENT + SCRNCOVID_bin,
            family = weibullPH)

#####################################################################################
# repeating everything for each prior mean used (log(2), log(1.2), log(1), log(3), log(4))
# at each mean used all phases are done and ppos is calculated

####### prior on HR  centered on log(2) #########
set.seed(123)
fit_HR2 <- fit.PH.HR(dataset = dat, 
                     brms_function1 = bbff1,
                     brms_function2 = bbff2,
                     HR.mean = log(2), 
                     HR.sd = sqrt(sd.fit)                  )

pred.datasets_HR2 <- predA(data_to_predict = dat0, 
                           data_observed = dat1, 
                           model1 = fit_HR2$Fit1, 
                           model2 = fit_HR2$Fit2, 
                           n.sim = nsim , cores)

fit.null_analysis <- compile.analysis.model(pred.datasets_HR2[[1]]). ## compiling model only the first time, then updating
analysis.datasets_HR2 <- mclapply(pred.datasets_HR2, \(x) main.analysis.update(x, fit.null_analysis), mc.cores = cores)
analysis.HR2 <- as.data.frame(do.call(rbind, analysis.datasets_HR2 %>% map(2)))

# ####### prior on HR  centered on log(1.2) #########
set.seed(123)
fit_HR1.2 <- fit.PH.HR(dataset = dat,
                       brms_function1 = bbff1,
                       brms_function2 = bbff2,
                       HR.mean = log(1.2),
                       HR.sd = sqrt(sd.fit)                  )

pred.datasets_HR1.2 <- predA(data_to_predict = dat0, 
                             data_observed = dat1, 
                             model1 = fit_HR1.2$Fit1, 
                             model2 = fit_HR1.2$Fit2, 
                             n.sim = nsim , cores)

analysis.datasets_HR1.2 <- mclapply(pred.datasets_HR1.2, \(x) main.analysis.update(x, fit.null_analysis), mc.cores = cores)
analysis.HR1.2 <- as.data.frame(do.call(rbind, analysis.datasets_HR1.2 %>% map(2)))

####### prior on HR  centered on log(1) #########
set.seed(123)
fit_HR1 <- fit.PH.HR(dataset = dat,
                     brms_function1 = bbff1,
                     brms_function2 = bbff2,
                     HR.mean = log(1),
                     HR.sd = sqrt(sd.fit)                  )

pred.datasets_HR1 <- predA(data_to_predict = dat0, 
                           data_observed = dat1, 
                           model1 = fit_HR1$Fit1, 
                           model2 = fit_HR1$Fit2, 
                           n.sim = nsim , cores)

analysis.datasets_HR1 <- mclapply(pred.datasets_HR1, \(x) main.analysis.update(x, fit.null_analysis), mc.cores = cores)
analysis.HR1 <- as.data.frame(do.call(rbind, analysis.datasets_HR1 %>% map(2)))

####### prior on HR  centered on log(3) #########
set.seed(123)
fit_HR3 <- fit.PH.HR(dataset = dat,
                     brms_function1 = bbff1,
                     brms_function2 = bbff2,
                     HR.mean = log(3),
                     HR.sd = sqrt(sd.fit)                  )

pred.datasets_HR3 <- predA(data_to_predict = dat0, 
                           data_observed = dat1, 
                           model1 = fit_HR3$Fit1, 
                           model2 = fit_HR3$Fit2, 
                           n.sim = nsim , cores)

analysis.datasets_HR3 <- mclapply(pred.datasets_HR3, \(x) main.analysis.update(x, fit.null_analysis), mc.cores = cores)
analysis.HR3 <- as.data.frame(do.call(rbind, analysis.datasets_HR3 %>% map(2)))

######## prior on HR  centered on log(4) #########
set.seed(123)
fit_HR4 <- fit.PH.HR(dataset = dat,
                     brms_function1 = bbff1,
                     brms_function2 = bbff2,
                     HR.mean = log(4),
                     HR.sd = sqrt(sd.fit))

pred.datasets_HR4 <- predA(data_to_predict = dat0, 
                           data_observed = dat1, 
                           model1 = fit_HR4$Fit1, 
                           model2 = fit_HR4$Fit2, 
                           n.sim = nsim , cores)

analysis.datasets_HR4 <- mclapply(pred.datasets_HR4, \(x) main.analysis.update(x, fit.null_analysis), mc.cores = cores)
analysis.HR4 <- as.data.frame(do.call(rbind, analysis.datasets_HR4 %>% map(2)))


##### saving ###
save.image(savedir)
