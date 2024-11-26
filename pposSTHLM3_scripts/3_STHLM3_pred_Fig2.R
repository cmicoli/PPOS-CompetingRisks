#==============================================================================
# FILENAME: 3_STHLM3_pred_Fig2.R
# PROJECT: 	STHLM3 PPOS
# PURPOSE: Predicting & analysis (PHASE 2 + 3 of PPOS) with modelling strategy A to produce Fig2 of manuscript
# AUTHOR: Chiara

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)
#==============================================================================

## PHASE 2/3 of PPOS: Predicting times after fitting models to STHLM3 data, then 
## censoring with 4 years of added follow-up and analysis the "final data", 
## saving 500 survivals to plot the results

########################################################################################
#  parameters that change for different analysis (number of sim/ cores/ FU time addes)
rm(list = ls())
print("model strategy A - diff splits")
load("/nfs/home/chimic/STHLM3/fit_pch_sthlm3_diffSplits_stan.RData")
splits <- "diff"
tchange1 <- c(0, tsplit1)
tchange2 <- c(0, tsplit2)
savedir <-  "/nfs/home/chimic/STHLM3/PRED_ToPLOT_sthlm3_diffSplits_stan_v2.RData" 
n.sim <- 2500
cores <- 40

source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
source("/nfs/home/chimic/STHLM3/00Functions_weibull_sim.R")
source("/nfs/home/chimic/STHLM3/00Functions_pch_sim.R")

########################################################################################
library(survival)
library(parallel)
library(rms)
library(MASS)
library(brms)
library(dplyr)
library(writexl)
library(rstanarm)
library(tidybayes)
library(tidyverse)
library(cmdstanr)
library(lubridate)
library(purrr)

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
                                          
                          predicted_cens <- cens_atT(predicted, 4 , predicted$maxFU18 ) ### point were you add time(extra years) to follow-up
                          predicted_dat_tot <- 
                            rbind(dat_inv1, dat_not.inv1) %>% mutate(event = NA, 
                                                                     time  = NA, 
                                                                     event.predicted = event18, 
                                                                     time.predicted = time18) %>%
                            #rbind(predicted)
                            rbind(predicted_cens)
                          
                          predicted_dat_tot

                        }, 
                        mc.cores = cores)

### Phase 3: Analysis - after censoring at different follow-up time
## analysis of 500 sampled predicted dataset to plot
predicted.surv_w <- mclapply(sample(1:n.sim, 500), \(x) surv.analysis.RR3(predicted_w[[x]]),
                             mc.cores = cores)

# saving results to plot
save( predicted.surv_w, file = savedir)



########################################################################################
#FIG2
person_main_surv$event18 <- as.factor(person_main_surv$event18)
aj.obs <- survfit(Surv(time18, event18) ~ group_rand, data = person_main_surv)
predicted.surv <- predicted.surv_w %>% map(1)
set.seed(123)

# Plot all the things adn saving it in "Output/survival_plot_241109.jpeg"
COL1 <- adjustcolor(c("grey", "transparent"), alpha.f = 0.10)
COL2 <- adjustcolor(c( "transparent", "grey"), alpha.f = 0.10)

jpeg("Output/survival_plot_241109.jpeg", width = 1400, height = 900, units = "px", quality = 300)  # Adjust dimensions and quality as needed
par(mfrow=c(2,2), mar=c(0.5, 0.5, 2, 0.2), oma = c(4, 5, 0.2, 0.2))
#1st row
plot(NULL, xlim = c(0, 11), ylim = c(0, 0.006), 
     pch=19, axes=F, type="b", cex=1.9, lwd=5,  
     xlab="Time (years)", ylab="Prostate Cancer mortality")
axis(2, cex.axis=1.3 , at = seq(0, 0.006, by = 0.001), labels = seq(0, 0.006, by = 0.001))
for (i in c(1:500)) {
  lines(predicted.surv[[i]], col = COL1, noplot = c("1", "(s0)"), lty = 1)
}
lines(aj.obs, conf.int = FALSE, lwd = 3, col = c("tomato4", "transparent"), noplot = c("1", "(s0)"), xaxs= "i", yaxs= "i")
title("Not invited group", cex.main = 1.9)
box(lty=1, col="black")

plot(NULL, xlim = c(0, 11), ylim = c(0, 0.006), xlab="Time (years)", ylab="Prostate Cancer mortality", 
     pch=19, axes=F, type="b", cex=1.9, lwd=5)
for (i in c(1:500)) {
  lines(predicted.surv[[i]], col = COL2, noplot = c("1", "(s0)"), lty = 1)
}
lines(aj.obs, conf.int = FALSE, lwd = 3, col = c("transparent", "navy"), noplot = c("1", "(s0)"), xaxs= "i", yaxs= "i")
title("Invited group", cex.main = 1.9)
box(lty=1, col="black")

# 2nd row
plot(NULL, xlim = c(0, 11), ylim = c(0, 0.12), xlab="Time (years)", ylab="Other-cause Mortality ", 
     pch=19, axes=F, type="b", cex=1.9, lwd=5)
axis(2, cex.axis=1.3 , at = seq(0, 0.12, by = 0.01), labels = seq(0, 0.12, by = 0.01))
axis(1, cex.axis=1.3, at = seq(0, 10, by = 1), labels = seq(0, 10, by = 1))
box(lty=1, col="black")
for (i in c(1:500)) {
  lines(predicted.surv[[i]], col = COL1, noplot = c("2", "(s0)"), lty = 1)
}
lines(aj.obs, conf.int = FALSE, lwd = 3, col = c("tomato4", "transparent"), noplot = c("2", "(s0)"), xaxs= "i", yaxs= "i")

plot(NULL, xlim = c(0, 11), ylim = c(0, 0.12), xlab="Time (years)", ylab="Other-cause Mortality", 
     pch=19, axes=F, type="b", cex.axis=1.9, lwd=5)
for (i in c(1:500)) {
  lines(predicted.surv[[i]], col = COL2, noplot = c("2", "(s0)"), lty = 1)
}
lines(aj.obs, conf.int = FALSE, lwd = 3, col = c("transparent", "navy"), noplot = c("2", "(s0)"), xaxs= "i", yaxs= "i")
axis(1, cex.axis=1.3, at = seq(0, 10, by = 1), labels = seq(0, 10, by = 1))
box(lty=1, col="black")

# Adding y-axis label for the first row
mtext("Prostate Cancer Mortality", side=2, line=2.5, cex=1.8, outer=TRUE, at=0.75)
# Adding y-axis label for the second row
mtext("Other Cause Mortality", side=2, line=2.5, cex=1.8, outer=TRUE, at=0.25)
# Adding x-axis label for the bottom row
mtext("Time since invitation (years)", side=1, line=2.5, cex=1.8, outer=TRUE, at=0.5)

# Save plot as JPEG
dev.off()

########################################################################################
sessionInfo()