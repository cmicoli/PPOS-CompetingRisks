
# Predictive probability of success (PPOS) - application with time-to-event data in the presence of competing event data. 

Chiara Micoli, Andrea Discacciati

Karolinska Institutet, Sweden

Published: 2024-08-12

Last updated: 2024-08-12

------------------------------------------------------------------------

Example for results presented in the project “Simulation-based Bayesian predictive probability of success for
interim monitoring of clinical trials with competing event data: two case studies”

------------------------------------------------------------------------

Functions can be `source`’d directly from R:

``` r
source("https://raw.githubusercontent.com/cmicoli/PPOS-CompetingRisks/main/00Functions_weibull_PPOS.R")
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
```

``` r
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

library(brms)
library(rstanarm)
library(survival)
library(tidybayes)
library(tidyverse)
library(parallel)
```

``` r
theme_set(theme_bw())
```

### Simulate right-censored survival data.

Simulating data until interim analysis, corresponding to an interim data
collection.

``` r
dat <- local({
  set.seed(1901)
  N <- 150
  TR <- rbinom(N, 1, 0.45)
  c <- rbinom(N, 1, 0.1)
  z <- rnorm(N)
  y <- flexsurv::rweibullPH(N, 
                            shape =  0.4, # gamma
                            scale = exp(-0.5 + log(0.5)*TR + log(0.2)*z)) # mu
  cens <- runif(N, 0, 60)
  time <- pmin(y, cens)
  status <- as.numeric(y <= cens)
  
  event <- ifelse(status == 1, rbinom(N, 1, 0.3) + 1, 0) #competing events 1-2
   
  data.frame(
    ID = 1:N,
    TIME = time,
    status = status,
    censored = 1 - status,
    EVENT = event,
    TREATMENT =  factor(TR,
                            levels = c(0, 1),   
                            labels = c("CNTR", "NEW_TR")), #factor(TR),
    c = factor(c),
    z = z
  )
})


head(dat)
```

     ID         TIME status censored EVENT TREATMENT c           z
      1 41.021730032      0        1     0      CNTR 0  1.37224254
      2  5.123673330      1        0     1      CNTR 0  0.24595045
      3 15.566988970      1        0     1      CNTR 0 -0.32276743
      4  0.006174074      1        0     2      CNTR 0  0.02744533
      5 27.295880378      0        1     0      CNTR 0  0.99959343
      6 42.601054255      0        1     0    NEW_TR 0  0.75863798

Presenting event summary for controls and new treatment group (85 CNTR
vs 65 NEW_TR): event 1 and event 2 are the two competing events, 0 for
censored data.

``` r
table(dat$TREATMENT,  useNA = "ifany")
```


      CNTR NEW_TR 
        85     65 

``` r
table(dat$TREATMENT, dat$EVENT, useNA = "ifany")
```

            
              0  1  2
      CNTR   22 46 17
      NEW_TR 32 22 11

``` r
summary(dat$TIME)
```

        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
     0.00001  0.05216  2.25706 10.64806 18.76344 59.80789 

Preparing data for fitting brms models.

``` r
dat$event.all <- ifelse(dat$EVENT %in% c(1,2), 1, 0)
dat$event1 <- dat$event.all * (dat$EVENT == 1)
dat$event2 <- dat$event.all * (dat$EVENT == 2)
dat$censored1 <- 1 - dat$event1
dat$censored2 <- 1 - dat$event2
```

## Definition of success

Time to recovery is analyzed in the trial by using a Weibull
distribution with prespecified weakly informative priors. To do so event
1 ($i=1$) and event 2 ($i=2$) are considered as competing events, and
two separate Bayesian proportional-hazard Weibull models will be used,
one for each competing outcome $i$. For each model, we will let the
logarithm of the Weibull distribution’s scale parameter u to depend on
treatment arm (X; new treatment or control arm) and on baseline value of
C, so that $log(u_i) = \alpha_i + \beta_i X + \gamma_i C$. The
distribution’s shape parameter ($v^i$), which is constrained to be
positive, will not depend on any covariates. The following priors are
used in the primary analysis:

$\alpha_i \sim Normal (0, 20)$

$\beta_i \sim Normal (0, \sqrt{0.5})$

$\gamma_i \sim Normal (0, \sqrt{0.5})$

$v_i \sim Exponential(1)$

Success is here defined as achieving the probability threshold for the
new treatment to graduate to the next phase of drug development. The
rule for graduation is formalized as
$P(csHR_{event1} > 1 | Data) \geq 0.975$.

## Analysis at interm

Start by fitting Fit Weibull PH model with brms. Model to fit: we
consider treatment arm and one binary covariate (c) in the models.

``` r
## models 
bbff1 <- bf(TIME | cens(censored1) ~   0 + Intercept + TREATMENT + c,
            family = weibullPH)
bbff2 <- bf(TIME | cens(censored2) ~   0 + Intercept + TREATMENT + c,
            family = weibullPH)
```

``` r
## output is median HR, and 95% qCrI for treatment effect, 
## and probability of recovery graduating : p(HR >1 | Dati) = # draws with HR >1 / 4000
fit.main.event1 <- brm(bbff1,
                  data = dat,     
                  chains = 4,
                  iter = 1500,
                  warmup = 500,
                  seed = 12345,
                  cores = 1,
                  stanvars = stanvars_weibullPH,
                  prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                            prior(normal(0, sqrt(.5)), class = b), 
                            prior(exponential(1), class = gamma)),
                  backend = "cmdstanr", ## rstan
                  refresh = 0)
```

    Running MCMC with 4 sequential chains...

    Chain 1 finished in 0.3 seconds.
    Chain 2 finished in 0.4 seconds.
    Chain 3 finished in 0.4 seconds.
    Chain 4 finished in 0.3 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.4 seconds.
    Total execution time: 2.0 seconds.

``` r
fit.main.event2 <- brm(bbff2,
                  data = dat,     
                  chains = 4,
                  iter = 1500,
                  warmup = 500,
                  seed = 12345,
                  cores = 1,
                  stanvars = stanvars_weibullPH,
                  prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                            prior(normal(0, sqrt(.5)), class = b), 
                            prior(exponential(1), class = gamma)),
                  backend = "cmdstanr", ## rstan
                  refresh = 0) 
```

    Running MCMC with 4 sequential chains...

    Chain 1 finished in 0.4 seconds.
    Chain 2 finished in 0.3 seconds.
    Chain 3 finished in 0.4 seconds.
    Chain 4 finished in 0.3 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.3 seconds.
    Total execution time: 1.7 seconds.

``` r
#fit.main.event1
pos.draws <- fit.main.event1 %>% 
    tidybayes::spread_draws(b_Intercept, b_TREATMENTNEW_TR , b_c1   , gamma)
  # print(dim(pos.draws))
  # print(colSums(is.na(pos.draws)))
HR.CrI <- quantile(exp(pos.draws$b_TREATMENTNEW_TR),
                     probs =c(0.5, 0.025, 0.975),  na.rm = FALSE) %>% 
    as.data.frame() %>% t() 
HR.CrI.pr <- HR.CrI %>% 
    cbind(count.HRvs1  = sum(exp(pos.draws$b_TREATMENTNEW_TR)>1)) %>% 
    as.data.frame() %>%
    mutate(pr = count.HRvs1 /dim(pos.draws)[1], ## probability of graduating in this dataset
           success = ifelse(pr > 0.975, 1, 0 )) ## definition of success with this dataset
HR.CrI.pr
```

            50%      2.5%     97.5% count.HRvs1     pr success
    . 0.5303063 0.3241296 0.8450462          10 0.0025       0

# Predictive probability of success analysis

### 1. Splitting observed data

Administrative censoring after 60 days.

``` r
# splitting between completely observed and censored data
dat1 <- dat %>% filter(event.all == 1 | TIME  == 60) %>% 
  select(ID, TREATMENT, c, TIME, EVENT)  ## 101 patients
dat0_cens <- dat %>% filter(event.all == 0 & TIME < 60) ## 32 patients
```

### 2. Defining number simulations K

``` r
nsim <- 3  ## K
# cores to use in the fitting if parallelized
cores <- 1 
```

### 3. Predicting baseline allocation of new patients

Maximum enrollement for the new treatment arm is 125.

``` r
extra.new.tr <- 125 - sum(dat$TREATMENT == "NEW_TR")
## new treatment has to be present 60 times, X is the number of controls needed to enroll 60 new_tr
D_new <- lapply (1:nsim, function(x) {
  n.control <- rnbinom(1, extra.new.tr, 0.45)
  d <- data.frame(
    ID = paste0("NEW", seq_len(n.control+ extra.new.tr)),
    TREATMENT = factor(rep(c("CNTR", "NEW_TR"), times = c(n.control, extra.new.tr))),
    c = factor(sample(c(0, 1), n.control+extra.new.tr, replace = TRUE, prob = c(0.9, 0.1)))
  )
})

# dat0 contains data that need prediction, either update of follow up for dat0_cens, or for D_new
dat0 <- lapply (D_new, function(x) {
  dat0 <- dat0_cens %>% 
    select(ID, TREATMENT, c, TIME, EVENT) %>%
    rbind (x %>% 
             mutate(TIME = 0, EVENT = 0))
})
```

Splitting into control and new treatment dataset, for future steps where
we will fit separate models

``` r
# new treatment datasets
dat_new.tr1 <- dat %>% filter(TREATMENT == "NEW_TR") %>% 
  filter(event.all == 1 | TIME  == 60) %>% 
  select(ID, TREATMENT, c, TIME, EVENT)  ## 45 patients
dat_new.tr0 <- lapply(dat0, \(x) x %>% filter(TREATMENT == "NEW_TR") )
                      #function(x) {  dat_new.tr0 <-  x %>% filter(TREATMENT == "NEW_TR")})
## 89 patients to predict: 32 already enrolled, 60 new

# control datasets
dat_control1 <- dat %>% filter(TREATMENT == "CNTR") %>% 
  filter(event.all == 1 | TIME  == 60) %>% 
  select(ID, TREATMENT, c, TIME, EVENT)  ## 55 patients
dat_control0 <- lapply(dat0, \(x) x %>% filter(TREATMENT == "CNTR") )
## XX patients to predict: 22 already enrolled, new patients in the trial (new accruals - from negative binomial)
```

## STRATEGY A

Modeling stratifying over treatment arm, using C as the only covariate
in the model.

### 4.A - Modeling step

``` r
## models 
bbff1.modelingA <- bf(TIME | cens(censored1) ~   0 + Intercept  + c,
            family = weibullPH)
bbff2.modelingA <- bf(TIME | cens(censored2) ~   0 + Intercept  + c,
            family = weibullPH)
```

``` r
# model fit for new treatment arm
fit_new.tr <- fit.noPH(dataset = dat %>% filter(TREATMENT == "NEW_TR"), 
                       brms_function1 = bbff1.modelingA,
                       brms_function2 = bbff2.modelingA   )
```

    Running MCMC with 4 parallel chains...

    Chain 1 finished in 0.4 seconds.
    Chain 2 finished in 0.4 seconds.
    Chain 3 finished in 0.4 seconds.
    Chain 4 finished in 0.4 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.4 seconds.
    Total execution time: 0.4 seconds.

    Running MCMC with 4 parallel chains...

    Chain 1 finished in 0.4 seconds.
    Chain 2 finished in 0.4 seconds.
    Chain 3 finished in 0.4 seconds.
    Chain 4 finished in 0.4 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.4 seconds.
    Total execution time: 0.5 seconds.

``` r
# model fit for control arm 
fit_control <- fit.noPH(dataset = dat %>% filter(TREATMENT == "CNTR"), 
                        brms_function1 = bbff1.modelingA,
                        brms_function2 = bbff2.modelingA   )
```

    Running MCMC with 4 parallel chains...

    Chain 1 finished in 0.5 seconds.
    Chain 2 finished in 0.5 seconds.
    Chain 3 finished in 0.5 seconds.
    Chain 4 finished in 0.5 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.5 seconds.
    Total execution time: 0.6 seconds.

    Running MCMC with 4 parallel chains...

    Chain 1 finished in 0.5 seconds.
    Chain 2 finished in 0.4 seconds.
    Chain 3 finished in 0.5 seconds.
    Chain 4 finished in 0.4 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.4 seconds.
    Total execution time: 0.6 seconds.

## 5.A - Predictive step

``` r
## new predicted dataset 
predicted_new.tr <- predicted_dataset_noPH(data_to_predict = dat_new.tr0, 
                                               data_observed = dat_new.tr1, 
                                               n.sim = nsim, 
                                               model1 = fit_new.tr$Fit1, 
                                               model2 = fit_new.tr$Fit2, 
                                               cores = cores)

predicted_control <- predicted_dataset_noPH(data_to_predict = dat_control0, 
                                                data_observed = dat_control1, 
                                                n.sim = nsim, 
                                                model1 = fit_control$Fit1, 
                                                model2 = fit_control$Fit2, 
                                                cores = cores)

#predicted <- lapply (1:nsim , \(x) rbind(predicted_control[[x]] , 
 #                                        predicted_new.tr [[x]]))

predicted <- lapply (1:nsim , \(x) rbind(predicted_control[[x]] , 
                                         predicted_new.tr [[x]]) %>%
                                    mutate(event.predicted = ifelse(is.na(event.predicted), 0,
                                                                    event.predicted), 
                                           time.predicted = ifelse(time.predicted == 0, 0.01,
                                                                   time.predicted)))
```

## 6.A - Analytic step

``` r
## analysis of new predicted dataset 
fit.null_analysis <- compile.analysis.model(predicted[[1]])
```

    Running MCMC with 1 chain...

    Chain 1 WARNING: No variance estimation is 
    Chain 1          performed for num_warmup < 20 
    Chain 1 Iteration: 1 / 2 [ 50%]  (Warmup) 
    Chain 1 Iteration: 2 / 2 [100%]  (Sampling) 
    Chain 1 finished in 0.0 seconds.

``` r
predicted_an <- mclapply(predicted, 
                         \(x)  main.analysis.update(dataset = x, fit.null = fit.null_analysis),
                         mc.cores = cores)
```

``` r
analysis.HR <- as.data.frame(do.call(rbind, predicted_an %>% map(2)))
head(analysis.HR)
```

       HR.median HR.025lower HR.975upper count.HRvs1      pr success
    .  0.4244328   0.2909498   0.6145362           0 0.00000       0
    .1 0.5209696   0.3668775   0.7410176           0 0.00000       0
    .2 0.5912076   0.4202998   0.8230841           5 0.00125       0

``` r
PPOS <- as.numeric(1 - (table(analysis.HR$success, useNA = "ifany") %>% prop.table() )[1])
PPOS
```

    [1] 0

## STRATEGY B

Modeling using PH assumption on the treatment effect and adding prior
belief on treatment effect, using C as extra covariate in the model.

Choosing $\mu$ and $\sigma$ for prior on the treatment effect

``` r
mu.fit <- 2
sd.fit <- 0.1 
```

### 4.B - Modeling step

Prior on treatment effect centered on $log(2)$

``` r
## models 
bbff1.modelingB <- bf(TIME | cens(censored1) ~   0 + Intercept + TREATMENT + c,
            family = weibullPH)
bbff2.modelingB <- bf(TIME | cens(censored2) ~   0 + Intercept + TREATMENT + c,
            family = weibullPH)
```

``` r
# model fit 
fit_HR2 <- fit.PH.HR(dataset = dat, 
                     brms_function1 = bbff1.modelingB,
                     brms_function2 = bbff2.modelingB,
                     HR.mean = log(mu.fit), 
                     HR.sd = sqrt(sd.fit))
```

    Running MCMC with 4 parallel chains...

    Chain 1 finished in 1.0 seconds.
    Chain 2 finished in 1.0 seconds.
    Chain 3 finished in 1.1 seconds.
    Chain 4 finished in 1.1 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 1.1 seconds.
    Total execution time: 1.1 seconds.

    Running MCMC with 4 parallel chains...

    Chain 2 finished in 0.8 seconds.
    Chain 1 finished in 0.9 seconds.
    Chain 3 finished in 0.9 seconds.
    Chain 4 finished in 0.8 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.9 seconds.
    Total execution time: 1.0 seconds.

### 5.B - Predictive step

``` r
pred.datasets_HR2 <- predicted_dataset_PH.HR(data_to_predict = dat0, 
                           data_observed = dat1, 
                           model1 = fit_HR2$Fit1, 
                           model2 = fit_HR2$Fit2, 
                           n.sim = nsim , cores)

pred.datasets_HR2 <- lapply (1:nsim , \(x) pred.datasets_HR2[[x]] %>%
                                    mutate(event.predicted = ifelse(is.na(event.predicted), 0,
                                                                    event.predicted), 
                                           time.predicted = ifelse(time.predicted == 0, 0.01, time.predicted)))
```

### 6.B - Analytic step

``` r
fit.null_analysis <- compile.analysis.model(pred.datasets_HR2[[1]])
```

    Running MCMC with 1 chain...

    Chain 1 WARNING: No variance estimation is 
    Chain 1          performed for num_warmup < 20 
    Chain 1 Iteration: 1 / 2 [ 50%]  (Warmup) 
    Chain 1 Iteration: 2 / 2 [100%]  (Sampling) 
    Chain 1 finished in 0.0 seconds.

``` r
analysis.datasets_HR2 <- mclapply(pred.datasets_HR2, 
                                  \(x) main.analysis.update(x, fit.null_analysis), mc.cores = cores)
```

``` r
analysis.HR2 <- as.data.frame(do.call(rbind, analysis.datasets_HR2 %>% map(2)))
```

``` r
table(analysis.HR2$success, useNA = "ifany")
```


    0 
    3 

``` r
PPOS <- as.numeric(1 - (table(analysis.HR2$success, useNA = "ifany") %>% prop.table() )[1])
PPOS
```

    [1] 0

## Session info.

``` r
sessionInfo()
```

    R version 4.3.2 (2023-10-31 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)

    Matrix products: default


    Random number generation:
     RNG:     L'Ecuyer-CMRG 
     Normal:  Inversion 
     Sample:  Rejection 
     
    locale:
    [1] LC_COLLATE=Swedish_Sweden.utf8  LC_CTYPE=Swedish_Sweden.utf8   
    [3] LC_MONETARY=Swedish_Sweden.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Swedish_Sweden.utf8    

    time zone: Europe/Stockholm
    tzcode source: internal

    attached base packages:
    [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
     [5] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
     [9] ggplot2_3.5.1   tidyverse_2.0.0 tidybayes_3.0.6 survival_3.5-8 
    [13] rstanarm_2.32.1 brms_2.21.0     Rcpp_1.0.12    

    loaded via a namespace (and not attached):
      [1] gridExtra_2.3        inline_0.3.19        rlang_1.1.3         
      [4] magrittr_2.0.3       matrixStats_1.3.0    compiler_4.3.2      
      [7] loo_2.7.0            vctrs_0.6.5          reshape2_1.4.4      
     [10] quadprog_1.5-8       pkgconfig_2.0.3      arrayhelpers_1.1-0  
     [13] fastmap_1.1.1        backports_1.4.1      utf8_1.2.4          
     [16] cmdstanr_0.7.1       threejs_0.3.3        promises_1.3.0      
     [19] deSolve_1.40         rmarkdown_2.26       markdown_1.12       
     [22] tzdb_0.4.0           ps_1.7.6             nloptr_2.0.3        
     [25] xfun_0.43            jsonlite_1.8.8       mstate_0.3.2        
     [28] later_1.3.2          R6_2.5.1             dygraphs_1.1.1.6    
     [31] stringi_1.8.3        StanHeaders_2.32.7   boot_1.3-28.1       
     [34] numDeriv_2016.8-1.1  rstan_2.32.6         knitr_1.46          
     [37] zoo_1.8-12           base64enc_0.1-3      bayesplot_1.11.1    
     [40] timechange_0.3.0     httpuv_1.6.15        Matrix_1.6-5        
     [43] splines_4.3.2        igraph_2.0.3         tidyselect_1.2.1    
     [46] rstudioapi_0.16.0    abind_1.4-5          yaml_2.3.8          
     [49] codetools_0.2-20     miniUI_0.1.1.1       processx_3.8.4      
     [52] curl_5.2.1           pkgbuild_1.4.4       lattice_0.22-6      
     [55] plyr_1.8.9           shiny_1.8.1.1        withr_3.0.0         
     [58] bridgesampling_1.1-2 posterior_1.5.0      coda_0.19-4.1       
     [61] evaluate_0.23        RcppParallel_5.1.7   muhaz_1.2.6.4       
     [64] ggdist_3.3.2         xts_0.13.2           pillar_1.9.0        
     [67] tensorA_0.36.2.1     checkmate_2.3.1      DT_0.33             
     [70] stats4_4.3.2         shinyjs_2.1.0        distributional_0.4.0
     [73] generics_0.1.3       hms_1.1.3            rstantools_2.4.0    
     [76] munsell_0.5.1        scales_1.3.0         minqa_1.2.6         
     [79] gtools_3.9.5         xtable_1.8-4         glue_1.7.0          
     [82] flexsurv_2.3         tools_4.3.2          shinystan_2.6.0     
     [85] data.table_1.15.4    lme4_1.1-35.3        colourpicker_1.3.0  
     [88] mvtnorm_1.2-4        grid_4.3.2           QuickJSR_1.1.3      
     [91] crosstalk_1.2.1      colorspace_2.1-0     nlme_3.1-164        
     [94] cli_3.6.2            fansi_1.0.6          svUnit_1.0.6        
     [97] Brobdingnag_1.2-9    V8_4.4.2             gtable_0.3.5        
    [100] digest_0.6.35        htmlwidgets_1.6.4    htmltools_0.5.8.1   
    [103] lifecycle_1.0.4      statmod_1.5.0        mime_0.12           
    [106] shinythemes_1.2.0    MASS_7.3-60.0.1     
