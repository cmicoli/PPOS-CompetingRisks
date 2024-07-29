# Predictive probability of success - application with time-to-event data in the presence of competing risks

Chiara Micoli

Karolinska Institutet, Sweden

Published: 2024-07-29

Last updated: 2024-07-29

------------------------------------------------------------------------

Example for results presented in the project *"Bayesian predictive probability of success for interim monitoring of clinical trials with competing risks data."*

Functions can be `source`’d directly from R:

``` r
source("https://raw.githubusercontent.com/anddis/brms-weibullPH/main/weibullPH_funs.R")
```

``` r
library(brms)
library(rstanarm)
library(survival)
library(tidybayes)
library(tidyverse)

theme_set(theme_bw())
```

### Simulate right-censored survival data.

``` r
simdata <- local({
  set.seed(1901)
  N <- 1000
  x <- rbinom(N, 1, 0.45)
  c <- rbinom(N, 1, 0.1)
  z <- rnorm(N)
  y <- flexsurv::rweibullPH(N, 
                            shape =  1.2, # gamma
                            scale = exp(0 + log(2)*x + log(0.75)*z)) # mu
  cens <- runif(N, 0, 4)
  time <- pmin(y, cens)
  status <- as.numeric(y <= cens)
  data.frame(
    time = time,
    status = status,
    censored = 1 - status,
    x = factor(x),
    c = factor(c),
    z = z
  )
})

head(simdata)
```

```
       time status censored x c           z
1 0.2560697      1        0 0 0 -1.72968590
2 1.7151737      1        0 0 0  0.89127671
3 0.6602768      1        0 0 0 -0.61354001
4 0.3558787      0        1 0 0 -0.41965989
5 0.5623408      1        0 0 0  0.02490977
6 0.1234761      1        0 1 0  0.34558304
```
