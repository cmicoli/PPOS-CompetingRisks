#==============================================================================
# FILENAME: 2a_STHLM3_pred&analysis_pch.R
# PROJECT: 	STHLM3 PPOS
# PURPOSE: Predicting & analysis (PHASE 2 + 3 of PPOS) with pch models fitted with STAN
# AUTHOR: Chiara Micoli, Andrea Discacciati

# R VERSION: R version 4.3.2 (2023-10-31)
            #Platform: x86_64-pc-linux-gnu (64-bit)
            #Running under: Red Hat Enterprise Linux 9.4 (Plow)
#==============================================================================

## PHASE 2/3 of PPOS: Predicting times after fitting models to STHLM3 data, then 
## censoring at different follow-up times and analysis the "final data", and recording PPOS value

########################################################################################
#  parameters that change for different analysis (number of sim/ cores/ FU time addes)
n.sim <- 2500
cores <- 40
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
source("/nfs/home/chimic/STHLM3/00Functions_weibull_sim.R")
source("/nfs/home/chimic/STHLM3/00Functions_pch_sim.R")
t.to.add <- c(4, 5, 6, 7, 8, 9, 10)  

########################################################################################
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
library(lubridate)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

########################################################################################
person_main_surv$event18 <- as.factor(person_main_surv$event18)
dat_tot <- person_main_surv[, c("LOPNR", "group_rand", "age_at_invitation", "date_invitation",
                                "time18", "event18", "event18.all" , "maxFU18"  , "Emigration_date"       )]
dat_tot$age_at_invitation_c <- dat_tot$age_at_invitation - mean(dat_tot$age_at_invitation)
summary(as.numeric(dat_tot$time18))

dat_inv <- dat_tot %>% filter(group_rand == "invited")
dat_not.inv <- dat_tot %>% filter(group_rand == "notinvited")


# splitting each dataset in 2 groups -> group1 don't need new prediction, group0 need to predict new times
dat_inv0 <- dat_inv  %>% 
  filter(event18.all == 0 & (is.na(Emigration_date) | Emigration_date> as.Date("2018-12-31")) )
dat_not.inv0 <- dat_not.inv %>% 
  filter(event18.all == 0 & (is.na(Emigration_date) | Emigration_date> as.Date("2018-12-31")) )


dat_inv1 <- dat_inv  %>% 
  filter(event18.all == 1 | Emigration_date <= as.Date("2018-12-31") )
dat_not.inv1 <- dat_not.inv %>% 
  filter(event18.all == 1 | Emigration_date <= as.Date("2018-12-31") )



# Phase 2: Prediction
# random beta draws depending on the model
if (length(tchange1) == 1 & length(tchange2) == 1) {
  r.beta_inv <- random.draw.beta.pch(fit1_inv, fit2_inv, n.sim)
  r.beta_not.inv <- random.draw.beta.pch(fit1_not.inv, fit2_not.inv, n.sim)
} else {
  r.beta_inv <- random.draw.beta.pch.stan(fit1_inv, fit2_inv, n.sim)
  r.beta_not.inv <- random.draw.beta.pch.stan(fit1_not.inv, fit2_not.inv, n.sim)
}
# Predicted datasets with times and event (censored at 6.5 y)
predicted_w <- mclapply(1:n.sim, 
                        function(x) {
                          predicted_inv <- rand_crisk_pch(
                            n = dim(dat_inv0)[1],
                            beta1 = as.numeric(r.beta_inv$fit.1[x,]),
                            beta2 = as.numeric(r.beta_inv$fit.2[x,]),
                            tchange1 , 
                            tchange2 ,
                            model_matrix = model.matrix(~ -1 + age_at_invitation_c, data = dat_inv0),
                            trunc = as.numeric(dat_inv0$time18) ### dat_inv0$maxFU18??? or time of real censoring
                          )
                          predicted_not.inv <- rand_crisk_pch(
                            n = dim(dat_not.inv0)[1],
                            beta1 = as.numeric(r.beta_not.inv$fit.1[x,]),
                            beta2 = as.numeric(r.beta_not.inv$fit.2[x,]),
                            tchange1 , 
                            tchange2 ,
                            model_matrix = model.matrix(~ -1 + age_at_invitation_c, data = dat_not.inv0),
                            trunc = as.numeric(dat_not.inv0$time18)
                          )
                          
                          predicted <- rbind(
                            cbind(dat_inv0 , predicted_inv), 
                            cbind(dat_not.inv0 , predicted_not.inv)
                          )
                          
                          predicted
                        }, 
                        mc.cores = cores)


### Phase 3: Analysis - after censoring at different follow-up time
## Results of prediction with different follow-up time added. 
# The results are saved in a list with the time added and the results of the prediction
# t.to.add <- c(4, 5, 6, 7, 8, 9, 10)  : has the possible amount of year to be added at the observed data
results_list <- list()
for (t in t.to.add) {
    timeadded <- t  # Update the timeadded variable for each iteration
# Predicted datasets with times and event (censored at maxFU18 + timeadded) -> to analyze with "surv.analysis.RR3"
    predicted.surv_w <- mclapply(predicted_w, 
                        function(x) {
                          predicted_cens <- cens_atT(x, timeadded , x$maxFU18 ) ### point were you add time(extra years) to follow-up
                          predicted_dat_tot <- 
                            rbind(dat_inv1, dat_not.inv1) %>% mutate(event = NA, 
                                                                     time  = NA, 
                                                                     event.predicted = event18, 
                                                                     time.predicted = time18) %>%
                            rbind(predicted_cens)
                          #predicted_dat_tot

                          surv.pred <- surv.analysis.RR3(predicted_dat_tot)
                          risk <- surv.pred[2:3]
                          risk
                        }, 
                        mc.cores = cores)

    risks_alpha_all <- as.data.frame(do.call(rbind, predicted.surv_w %>% map(1)))

    ppos3.5 <- sum(as.numeric(risks_alpha_all$p.value_PCSM )< 0.035)/length(risks_alpha_all$p.value_PCSM )
    
    # Save the results for the current t in the results_list
    results_list[[t]] <- list(ppos0.035 = ppos3.5, risks_alpha_all = risks_alpha_all)                           
}

# Create a data frame with the results
ppos_t <- cbind(t.to.add, end_FU_time = as.Date("2018-12-31") + years(t.to.add), 
        as.data.frame(do.call(rbind, results_list %>% map(1)))) %>% 
  rename("time" = "t.to.add") 

ppos_t

#saving results
save( results_list, ppos_t, file=savedir)
