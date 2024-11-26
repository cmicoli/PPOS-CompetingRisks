#==============================================================================
# FILENAME: 3_sensitivity_plot.R
# PROJECT: 	ISPY-COVID PPOS
# PURPOSE: Plot of the Sensitivity PPOS analysis  - FIG1 of the manuscript
# AUTHOR: Chiara Micoli, Andrea Discacciati

# R VERSION: R version 4.3.2 (2023-10-31)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Red Hat Enterprise Linux 9.4 (Plow)
#==============================================================================
## Plotting the sensitivity analysis for the prior on the treatment effect in the weibullPH models

library(survival)
library(parallel)
library(rms)
library(MASS)
library(brms)
library(bbmle)
library(dplyr)
library(huxtable)
library(purrr)
library(forcats)
library(ggplot2)
library(latex2exp)
library(tidyverse)
################################################################################
# functions that will be needed to merge the data
summary.ppos <- function (colname) {
  df <- lapply(list(analysis.HR1,   analysis.HR1.2, 
                    analysis.HR2, analysis.HR3, analysis.HR4 ),
               function(HR) {
                 as.numeric((table(HR$success, useNA = "ifany") %>% prop.table() )[2])}) 
  df <- as.data.frame(do.call(rbind, df))
  rownames(df) <- c("HR1", "HR1.2", "HR2", "HR3", "HR4")
  colnames(df) <- c(colname)
  return(df)}

summary.ppos2 <- function (sd) {
  df <- lapply(list(analysis.HR1,   analysis.HR1.2, 
                    analysis.HR2, analysis.HR3, analysis.HR4 ),
               function(HR) {
                 as.numeric((table(HR$success, useNA = "ifany") %>% prop.table() )[2])}) 
  df <- as.data.frame(do.call(rbind, df))
  df$prior.HR <- c(1, 1.2, 2, 3, 4)
  df$prior.SD <- sd
  rownames(df) <- c("HR1", "HR1.2", "HR2", "HR3", "HR4")
  colnames(df) <- c("ppos", "prior.HR", "prior.SD")
  return(df)}

join.dat <- function (sd , 
                      analysis.HR1 = analysis.HR1 ,   
                      analysis.HR1.2 = analysis.HR1.2,
                      analysis.HR2 = analysis.HR2, 
                      analysis.HR3 = analysis.HR3, 
                      analysis.HR4 = analysis.HR4 ) {
  
  analysis.HR1$prior.HR <- 1
  analysis.HR1.2$prior.HR <- 1.2
  analysis.HR2$prior.HR <- 2
  analysis.HR3$prior.HR <- 3
  analysis.HR4$prior.HR <- 4
  
  analysis.HR.tot <- rbind(analysis.HR1,   analysis.HR1.2, 
                           analysis.HR2, analysis.HR3, analysis.HR4 ) %>% 
    mutate(prior.SD = sd)
  return(analysis.HR.tot)
}
################################################################################
# loading different results
load("Output/vector/fit_predictions_priorPH_sd05_A2500_20241106.RData")
PPOS.05 <- summary.ppos("sd05")
PPOS.05.long <- summary.ppos2(0.5)

dat.05 <- join.dat(sd = 0.5, 
                   analysis.HR1,   analysis.HR1.2, 
                   analysis.HR2, analysis.HR3, analysis.HR4 )

load("Output/vector/fit_predictions_priorPH_sd02_A2500_20241106.RData")
PPOS.02 <- summary.ppos("sd02")
PPOS.02.long <- summary.ppos2(0.2)

dat.02 <- join.dat(sd = 0.2, 
                   analysis.HR1,   analysis.HR1.2, 
                   analysis.HR2, analysis.HR3, analysis.HR4 )

load("Output/vector/fit_predictions_priorPH_sd01_A2500_20241106.RData")
PPOS.01 <- summary.ppos("sd01")
PPOS.01.long <- summary.ppos2(0.1)

dat.01 <- join.dat(sd = 0.1, 
                   analysis.HR1,   analysis.HR1.2, 
                   analysis.HR2, analysis.HR3, analysis.HR4 )

################################################################################
# formatting the data to be then plotted
PPOS <- cbind (PPOS.05, PPOS.02, PPOS.01)
PPOS
PPOS.long <- rbind (PPOS.05.long, PPOS.02.long, PPOS.01.long) %>% 
  mutate(ppos.perc = as.numeric(num(ppos*100 , digits = 2)), 
         ppos = round(ppos , digits = 3))
PPOS.long

PPOS.long$label <- sprintf(
  "\U03BC = log(%s), \U03C3 = %s",
  PPOS.long$prior.HR,
  PPOS.long$prior.SD)
PPOS.long$label.p <- sprintf(
  "%.2f%%", ## "PPOS = %s %%"
  PPOS.long$ppos.perc)

PPOS.long$p.HR <- sprintf(
  "\U03BC = log(%s)",
  PPOS.long$prior.HR)
PPOS.long$p.SD <- sprintf(
  "\U03C3 = √%s",
  PPOS.long$prior.SD)


analysis.HR.tot <- rbind(  dat.05, dat.02, dat.01)
analysis.HR.tot$p.HR <- sprintf(
  "\U03BC = log(%s)",  ## mu = U03BC
  analysis.HR.tot$prior.HR)
analysis.HR.tot$p.SD <- sprintf(
  "\U03C3 = √%s",  ## sigma  = U03C3
  analysis.HR.tot$prior.SD)

analysis.HR.tot$prior.SD <- fct_rev(factor(analysis.HR.tot$prior.SD ))
analysis.HR.tot$p.SD <- factor(analysis.HR.tot$p.SD , levels = c("\U03C3 = √0.5", "\U03C3 = √0.2", "\U03C3 = √0.1"))
PPOS.long$p.SD <- fct_rev(factor(PPOS.long$p.SD ))

################################################################################
plot.fig1 <- ggplot(data = analysis.HR.tot , ##%>% filter(prior.HR != 4) ,
                       aes(x=pr)) + 
  geom_histogram(aes( y=after_stat(density*width), 
                      fill = after_stat(x>0.975), 
                      colour = after_stat(x>0.975)   ## after_stat(cut(x, breaks))
  ), 
  lwd = 0.25, breaks = seq(0, 1, by = 0.025), alpha=0.5)+   #binwidth= 0.05, 
  labs(x = expression("Pr(" ~ "csHR"["recovery"] ~ ">1 |" ~ tilde("D")["f"] ~ ")"),   #"Pr(HR recovery >1 | data)",
       y = "Proportion") +
  scale_y_continuous( 
    expand = c(0,0.01) 
  ) + 
  scale_x_continuous(expand = c(0.01,0), labels = prettyNum(c(0,0.25,0.50, 0.75,1)) ) + 
  scale_color_manual(values=c("black", "tomato4"))+
  scale_fill_manual(values=c("grey", "tomato4"))+
  theme_bw()+theme(legend.position = "none", 
                   text=element_text(size=16), 
                   panel.grid.major = element_blank(),  # Remove major gridlines
                   panel.grid.minor = element_blank()) +
  facet_grid(p.HR ~ p.SD) +   #p.SD ~ p.HR , switch = "y"
  geom_text(
    size    = 5, 
    data    = PPOS.long,
    mapping = aes(x = 0.97, y = 0.5, label = sprintf(ppos, fmt = '%#.3f')),
    hjust   = 1.05,
    vjust   = 1.5, color = "tomato4"
  )
## remove switch and position in scale_y_continous to have axis on the left side
plot.fig1