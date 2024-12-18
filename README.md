
# Predictive probability of success (PPoS) - application to time-to-event data in the presence of competing risks

Chiara Micoli, Andrea Discacciati

Karolinska Institutet, Sweden

Published: 2024-08-12

Last updated: 2024-12-18

------------------------------------------------------------------------

Reproducible example using simulated data for the the paper
“Simulation-based Bayesian predictive probability of success for interim
monitoring of clinical trials with competing event data: two case
studies”.

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
library(survival)
library(flexsurv)
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
  c <- rbinom(N, 1, 0.2)
  z <- rnorm(N)
  y <- rweibullPH(N, 
                  shape =  1.5, # gamma
                  scale = exp(-6 + log(1.4)*TR + log(0.8)*z)) # mu
  cens <- runif(N, 0, 60)
  time <- pmin(y, cens)
  status <- as.numeric(y <= cens)
  
  event <- ifelse(status == 1, rbinom(N, 1, 0.3)+1, 0) #competing events 1-2
  
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

      ID      TIME status censored EVENT TREATMENT c           z
    1  1 41.021730      0        1     0      CNTR 1  1.37224254
    2  2 47.274011      0        1     0      CNTR 0  0.24595045
    3  3 37.410411      0        1     0      CNTR 0 -0.32276743
    4  4  9.822373      1        0     2      CNTR 1  0.02744533
    5  5 27.295880      0        1     0      CNTR 1  0.99959343
    6  6 42.601054      0        1     0    NEW_TR 0  0.75863798

Presenting event summary for controls and new treatment group: event 1
and event 2 are the two competing events, 0 for censored data.

``` r
table(dat$TREATMENT,  useNA = "ifany")
```


      CNTR NEW_TR 
        85     65 

``` r
table(dat$TREATMENT, dat$EVENT, useNA = "ifany")
```

            
              0  1  2
      CNTR   51 21 13
      NEW_TR 43 17  5

``` r
summary(dat$TIME)
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     0.5683 12.8427 20.8761 22.8823 32.9039 54.3467 

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

    Chain 1 finished in 0.6 seconds.
    Chain 2 finished in 0.7 seconds.
    Chain 3 finished in 0.6 seconds.
    Chain 4 finished in 0.6 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.6 seconds.
    Total execution time: 2.7 seconds.

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

    Chain 1 finished in 0.7 seconds.
    Chain 2 finished in 0.6 seconds.
    Chain 3 finished in 0.5 seconds.
    Chain 4 finished in 0.5 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.6 seconds.
    Total execution time: 2.5 seconds.

``` r
#fit.main.event1
pos.draws <- fit.main.event1 %>% 
  spread_draws(b_Intercept, b_TREATMENTNEW_TR , b_c1   , gamma)
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

           50%      2.5%    97.5% count.HRvs1    pr success
    . 1.265445 0.7109435 2.281368        3156 0.789       0

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
nsim <- 100  ## K
# cores to use in the fitting if parallelized
cores <- 4
```

### 3. Predicting baseline allocation of new patients

Maximum enrolment for the new treatment arm is 125.

``` r
extra.new.tr <- 125 - sum(dat$TREATMENT == "NEW_TR")
## new treatment has to be present 60 times, X is the number of controls needed to enroll 60 new_tr
D_new <- lapply (1:nsim, function(x) {
  n.control <- rnbinom(1, extra.new.tr, 0.45)
  d <- data.frame(
    ID = paste0("NEW", seq_len(n.control+ extra.new.tr)),
    TREATMENT = factor(rep(c("CNTR", "NEW_TR"), times = c(n.control, extra.new.tr))),
    c = factor(rbinom(n.control+extra.new.tr, 1, 
                      rbeta(1, table(dat$c)["1"]+1, sum(table(dat$c))-table(dat$c)["1"]+1)),
               levels = c(0, 1))
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

## Main analysis (A)

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

    Chain 1 finished in 0.7 seconds.
    Chain 2 finished in 0.6 seconds.
    Chain 3 finished in 0.6 seconds.
    Chain 4 finished in 0.6 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.6 seconds.
    Total execution time: 0.8 seconds.

    Running MCMC with 4 parallel chains...

    Chain 1 finished in 0.6 seconds.
    Chain 2 finished in 0.5 seconds.
    Chain 3 finished in 0.5 seconds.
    Chain 4 finished in 0.5 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.5 seconds.
    Total execution time: 0.7 seconds.

``` r
# model fit for control arm 
fit_control <- fit.noPH(dataset = dat %>% filter(TREATMENT == "CNTR"), 
                        brms_function1 = bbff1.modelingA,
                        brms_function2 = bbff2.modelingA   )
```

    Running MCMC with 4 parallel chains...

    Chain 2 finished in 0.7 seconds.
    Chain 1 finished in 0.9 seconds.
    Chain 3 finished in 0.8 seconds.
    Chain 4 finished in 0.8 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.8 seconds.
    Total execution time: 1.0 seconds.

    Running MCMC with 4 parallel chains...

    Chain 2 finished in 0.7 seconds.
    Chain 4 finished in 0.7 seconds.
    Chain 1 finished in 0.9 seconds.
    Chain 3 finished in 0.8 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 0.8 seconds.
    Total execution time: 1.0 seconds.

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

predicted <- lapply (1:nsim , \(x) rbind(predicted_control[[x]] , 
                                         predicted_new.tr[[x]]))
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
analysis.HR$success <- factor(analysis.HR$success, levels = 0:1)
head(analysis.HR)
```

       HR.median HR.025lower HR.975upper count.HRvs1      pr success
    .  1.0289299   0.7478691    1.423476        2287 0.57175       0
    .1 0.7960582   0.5766778    1.110761         352 0.08800       0
    .2 1.2763156   0.9460287    1.693957        3782 0.94550       0
    .3 1.0912287   0.8143284    1.462149        2847 0.71175       0
    .4 1.8410029   1.3340284    2.498667        4000 1.00000       1
    .5 1.3031306   0.9635078    1.767451        3820 0.95500       0

``` r
PPOS <- (table(analysis.HR$success, useNA = "ifany") %>% prop.table())["1"]
PPOS
```

       1 
    0.44 

## Sensitivity analysis (B)

Modeling using PH assumption on the treatment effect and adding prior
belief on treatment effect, using C as extra covariate in the model.

Choosing $\mu$ (mu.fit) and $\sigma^2$ (var.fit) for prior on the treatment effect

``` r
mu.fit <- 2
var.fit <- 0.1 
```

### 4.B - Modeling step

Prior on log-treatment effect (log(csHR)) centered on $log(2)$

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
                     logHR.mean = log(mu.fit), 
                     logHR.sd = sqrt(var.fit))
```

    Running MCMC with 4 parallel chains...

    Chain 2 finished in 1.4 seconds.
    Chain 3 finished in 1.4 seconds.
    Chain 1 finished in 1.4 seconds.
    Chain 4 finished in 1.4 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 1.4 seconds.
    Total execution time: 1.5 seconds.

    Running MCMC with 4 parallel chains...

    Chain 2 finished in 1.3 seconds.
    Chain 3 finished in 1.3 seconds.
    Chain 4 finished in 1.2 seconds.
    Chain 1 finished in 1.5 seconds.

    All 4 chains finished successfully.
    Mean chain execution time: 1.3 seconds.
    Total execution time: 1.6 seconds.

### 5.B - Predictive step

``` r
pred.datasets_HR2 <- predicted_dataset_PH.HR(data_to_predict = dat0, 
                                             data_observed = dat1, 
                                             model1 = fit_HR2$Fit1, 
                                             model2 = fit_HR2$Fit2, 
                                             n.sim = nsim , cores)
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
analysis.HR2$success <- factor(analysis.HR2$success, levels = 0:1)
```

``` r
table(analysis.HR2$success, useNA = "ifany")
```


     0  1 
    25 75 

``` r
PPOS <- (table(analysis.HR2$success, useNA = "ifany") %>% prop.table())["1"]
PPOS
```

       1 
    0.75 

# Session info

``` r
sessionInfo()
```

    R version 4.4.2 (2024-10-31)
    Platform: aarch64-apple-darwin20
    Running under: macOS Sequoia 15.1

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

    Random number generation:
     RNG:     L'Ecuyer-CMRG 
     Normal:  Inversion 
     Sample:  Rejection 
     
    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: Europe/Stockholm
    tzcode source: internal

    attached base packages:
    [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] rstan_2.32.6        StanHeaders_2.32.10 lubridate_1.9.3    
     [4] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
     [7] purrr_1.0.2         readr_2.1.5         tidyr_1.3.1        
    [10] tibble_3.2.1        ggplot2_3.5.1       tidyverse_2.0.0    
    [13] tidybayes_3.0.7     flexsurv_2.3.2      survival_3.7-0     
    [16] brms_2.22.0         Rcpp_1.0.13-1      

    loaded via a namespace (and not attached):
     [1] tidyselect_1.2.1     svUnit_1.0.6         loo_2.8.0           
     [4] fastmap_1.2.0        TH.data_1.1-2        tensorA_0.36.2.1    
     [7] digest_0.6.37        estimability_1.5.1   timechange_0.3.0    
    [10] lifecycle_1.0.4      statmod_1.5.0        processx_3.8.4      
    [13] magrittr_2.0.3       posterior_1.6.0      compiler_4.4.2      
    [16] rlang_1.1.4          tools_4.4.2          utf8_1.2.4          
    [19] yaml_2.3.10          data.table_1.16.2    knitr_1.49          
    [22] bridgesampling_1.1-2 pkgbuild_1.4.5       cmdstanr_0.8.0      
    [25] multcomp_1.4-26      abind_1.4-8          withr_3.0.2         
    [28] numDeriv_2016.8-1.1  stats4_4.4.2         grid_4.4.2          
    [31] mstate_0.3.3         fansi_1.0.6          inline_0.3.20       
    [34] xtable_1.8-4         colorspace_2.1-1     emmeans_1.10.5      
    [37] scales_1.3.0         MASS_7.3-61          cli_3.6.3           
    [40] mvtnorm_1.3-2        rmarkdown_2.29       generics_0.1.3      
    [43] RcppParallel_5.1.9   rstudioapi_0.17.1    tzdb_0.4.0          
    [46] splines_4.4.2        bayesplot_1.11.1     muhaz_1.2.6.4       
    [49] matrixStats_1.4.1    vctrs_0.6.5          Matrix_1.7-1        
    [52] sandwich_3.1-1       jsonlite_1.8.9       hms_1.1.3           
    [55] arrayhelpers_1.1-0   ggdist_3.3.2         glue_1.8.0          
    [58] codetools_0.2-20     ps_1.8.1             distributional_0.5.0
    [61] stringi_1.8.4        gtable_0.3.6         QuickJSR_1.4.0      
    [64] quadprog_1.5-8       munsell_0.5.1        pillar_1.9.0        
    [67] htmltools_0.5.8.1    Brobdingnag_1.2-9    deSolve_1.40        
    [70] R6_2.5.1             evaluate_1.0.1       lattice_0.22-6      
    [73] backports_1.5.0      rstantools_2.4.0     gridExtra_2.3       
    [76] coda_0.19-4.1        nlme_3.1-166         checkmate_2.3.2     
    [79] xfun_0.49            zoo_1.8-12           pkgconfig_2.0.3     
