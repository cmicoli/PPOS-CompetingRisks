## functions for pch models  
#
# [hHpdq]pch functions are from PWEALL v1.3.0.1 (CRAN). Translated from Fortran
# to R.
#
# beta1 and beta2 contain log-baseline hazards (and log(HRs)). If a dummy
# parametrisation without intercept for the PCH model was used, then no
# pre-processing to the parameter vectors is necessary. Otherwise things like
# coef(pch)[1]+c(0, coef(pch)[2:x]) for standard contrast or
# cumsum(coef(pch)[1:x]) for sequential difference parametrisation are
# necessary.
#
# brms exponential distribution uses a 1/lambda=exp(xb) parametrisation,
# coefficients' vector must be multiplied by -1.

random.draw.beta.pch.stan <- function(fit1, fit2, n.sim){
  b1 <- fit1$draws(format = "data.frame")
  b2 <- fit2$draws(format = "data.frame")
  fit <- list(
    fit.1 = - b1[, grepl("^b", names(b1))],  ## for pch models in brm (stan) invert sign of parameters
    fit.2 = - b2[, grepl("^b", names(b2))]
  )
  
  r.beta <- lapply(fit, function(x) {x[do.call("sample", list(nrow(x),n.sim)),]})  
 # r.beta <- lapply(r.beta, function(x) {-as.numeric(x)})
  
  r.beta
}

# to use if models are fitted with brms
random.draw.beta.pch <- function(fit1, fit2, n.sim) {
  b1 <- as.data.frame(fit1)
  b2 <- as.data.frame(fit2)
  fit <- list(
    fit.1 = - (b1[, grepl("^b", names(b1))]),  ## for pch models in brm (stan) invert sign of parameters
    fit.2 = - (b2[, grepl("^b", names(b2))])
  )
  
  r.beta <- lapply(fit, function(x) {x[do.call("sample", list(nrow(x),n.sim)),]})  
  #r.beta <- lapply(r.beta, function(x) {-as.numeric(x)})
  r.beta
}


## hazard within time-window
hpch <- function(t, hazard, tchange) {
  stopifnot(t >= 0)
  stopifnot(tchange[1] == 0)
  stopifnot(length(hazard) == length(tchange))
  stopifnot(hazard >= 0)
  #   PWEALL::pwe(t = t, rate = hazard, tchange = tchange)$hazard
  
  hazard[findInterval(t, tchange)]
}



## cumulative hazard within time-window
Hpch <- function(t, hazard, tchange) {
  stopifnot(t >= 0)
  stopifnot(tchange[1] == 0)
  stopifnot(length(hazard) == length(tchange))
  stopifnot(hazard >= 0)
  # PWEALL::pwe(t = t, rate = hazard, tchange = tchange)$cumhazard
  
  i <- findInterval(t, tchange)
  
  vapply(
    seq_len(length(t)), function(j) {
      
      if (i[j] == 1) {
        hazard[i[j]] * t[j]
      } else {
        length.hazard <- length(hazard)
        tmax <- max(t[j], max(tchange))+1
        tplus <- rep(0, length.hazard)
        tplus[length.hazard] <- tmax
        tplus[1:(length.hazard-1)] <- tchange[2:length.hazard]
        rt <- hazard * (tplus - tchange)
        
        sum(rt[1:(i[j]-1)])+hazard[i[j]]*(t[j] - tchange[i[j]])
      } # end if
    }, # end function
    FUN.VALUE = numeric(1)
  )
}
ppch <- function(t, hazard, tchange, lower.tail = TRUE) {
  # PWEALL::pwe(t = t, rate = hazard, tchange = tchange)$surv
  
  S <- exp(-Hpch(t, hazard, tchange))
  if (lower.tail)
    return(1-S)
  return(S)
}


dpch <- function(t, hazard, tchange) {
  # PWEALL::pwe(t = t, rate = hazard, tchange = tchange)$density
  
  ppch(t, hazard, tchange, lower.tail = FALSE) * hpch(t, hazard, tchange)
}


qpch <- function (p, hazard, tchange) {
  # PWEALL::qpwe(p = p, rate = hazard, tchange = tchange
  np <- length(p)
  nr <- length(hazard)
  
  ss <- rep(0, nr)
  if (nr >= 2) 
    for(i in seq_len(nr-1)) 
      ss[i+1] = ss[i] + hazard[i]*(tchange[i+1]-tchange[i])
  
  ap <- -log1p(-p)
  
  outr <- vapply(seq_len(np), \(i) {
    j <- sum(ss <= ap[i])
    if (j == 0)
      0
    else
      tchange[j] + (ap[i]-ss[j])/hazard[j]
  }, FUN.VALUE = numeric(1))
  
  outr
}


quantilef_all_pch <- function (beta1, 
                               beta2, 
                               tchange1, 
                               tchange2, 
                               model_matrix, 
                               trunc) {
  stopifnot(length(trunc) == nrow(model_matrix))
  
  n.pred <- nrow(model_matrix)
  
  n.pieces1 <- length(tchange1)
  n.pieces2 <- length(tchange2)
  
  length.beta1 <- length(beta1)
  length.beta2 <- length(beta2)
  
  tchange.union <- sort(base::union(tchange1, tchange2))
  
  if (isTRUE(tchange.union == 0)) { # constant hazard
  hazard1 <- exp(model_matrix %*% beta1[(n.pieces1+1):length.beta1]) * exp(beta1[1])
  hazard2 <- exp(model_matrix %*% beta2[(n.pieces2+1):length.beta2]) * exp(beta2[1])
  hazard.sum <- hazard1 + hazard2
  
  cdf.at.trunc <- ifelse(trunc == 0, 0, pexp(q = trunc, rate = hazard.sum))
  u <- runif(n.pred, min = cdf.at.trunc, max = 1)
  t <- qexp(p = u, rate = hazard.sum)

} else { # piecewise constant hazard
  hazard1 <- exp(model_matrix %*% beta1[(n.pieces1+1):length.beta1]) %*% hpch(t = tchange1, hazard = exp(beta1[1:n.pieces1]), tchange = tchange1)
  hazard2 <- exp(model_matrix %*% beta2[(n.pieces2+1):length.beta2]) %*% hpch(t = tchange2, hazard = exp(beta2[1:n.pieces2]), tchange = tchange2)
  hazard.sum    <- t(apply(hazard1, 1, \(x) hpch(tchange.union, x, tchange1)) + apply(hazard2, 1, \(x) hpch(tchange.union, x, tchange2)))
  
  t <- vapply(seq_len(n.pred), \(j) {
    cdf.at.trunc <- if (trunc[j] == 0) 0 else ppch(trunc[j], hazard = hazard.sum[j, ], tchange = tchange.union)
    u <- runif(1, min = cdf.at.trunc, max = 1)
    qpch(p = u, hazard = hazard.sum[j, ], tchange = tchange.union)
  }, FUN.VALUE = numeric(1))
}
  # OBS: hazard[12] = HR * baseline_hazard
  #hazard1 <- exp(model_matrix %*% beta1[(n.pieces1+1):length.beta1]) %*% hpch(t = tchange1, hazard = exp(beta1[1:n.pieces1]), tchange = tchange1)
  #hazard2 <- exp(model_matrix %*% beta2[(n.pieces2+1):length.beta2]) %*% hpch(t = tchange2, hazard = exp(beta2[1:n.pieces2]), tchange = tchange2)
  # hazard.sum    <- t(apply(hazard1, 1, \(x) hpch(tchange.union, x, tchange1)) + apply(hazard2, 1, \(x) hpch(tchange.union, x, tchange2)))
  #if (n.pieces1 > 1 | n.pieces2 > 1) {
   # hazard.sum <- t(apply(hazard1, 1, \(x) hpch(tchange.union, x, tchange1)) + apply(hazard2, 1, \(x) hpch(tchange.union, x, tchange2)))
  #} else {
  #  hazard.sum <- (as.matrix(apply(hazard1, 1, \(x) hpch(tchange.union, x, tchange1))) + as.matrix(apply(hazard2, 1, \(x) hpch(tchange.union, x, tchange2))))
  #}
  
  
  #t <- vapply(seq_len(n.pred), \(j) {
   # cdf.at.trunc <- if (trunc[j] == 0) 0 else ppch(trunc[j], hazard = hazard.sum[j, ], tchange = tchange.union)
    
   # u <- runif(1, min = cdf.at.trunc, max = 1)
   # qpch(p = u, hazard = hazard.sum[j, ], tchange = tchange.union)
  #}, FUN.VALUE = numeric(1))
  
  return(t)
}

rand_crisk_pch <- function(n, beta1, beta2, tchange1, tchange2, model_matrix, trunc) {
  n.pieces1 <- length(tchange1)
  n.pieces2 <- length(tchange2)
  
  length.beta1 <- length(beta1)
  length.beta2 <- length(beta2)
  
  # step 2 of Beyersmann (from u to random event times - we don't know what competing event happened yet)
  event_times <- quantilef_all_pch(beta1 = beta1, beta2 = beta2, tchange1 = tchange1, tchange2 = tchange2, model_matrix = model_matrix, trunc)
  
  # step 3 of Beyersmann (which competing event happened at each event time?)
  hazard1 <- hpch(t = event_times, hazard = exp(beta1[1:n.pieces1]), tchange = tchange1) * exp(model_matrix %*% beta1[(n.pieces1+1):length.beta1])
  hazard2 <- hpch(t = event_times, hazard = exp(beta2[1:n.pieces2]), tchange = tchange2) * exp(model_matrix %*% beta2[(n.pieces2+1):length.beta2])
  
  causes <- rbinom(n, 1, prob = hazard2/(hazard1+hazard2)) + 1  
  
  data.frame(
    time = event_times,
    event = causes  
  )
}




## function data add censorship at time t
cens_atT <- function(predicted,## predicted dataset
                     t_added, ## t_added = how much follow-up time do I want to simulate
                     maxFU18) { # max follow-up at censoring 
  #dat){ ## dataset of censored obs that were used in rand_crisk 
  predicted$time.predicted <- ifelse(predicted$time <= maxFU18+t_added, predicted$time, maxFU18 + t_added)
  predicted$event.predicted <- factor(ifelse(predicted$time <= maxFU18 + t_added, predicted$event, "0"))
  
  predicted
}
# predicted output of this function is a data set with censored predicted time




## function that from RR conf interval tells if is significant or not (1 = risk reduction, 0 = no difference, 2 = higher risk)
alpha_sign_v2 <- function (rr){
  # 1 for significant  reduction result, 0 for non significant result or higher risk for invited
  alpha <- sapply(1:length(rr$cause), function(x){ifelse(rr[x,"CI95.1"]< 1 & rr[x, "CI95.2"]< 1, 1, 
                                                         ifelse(rr[x,"CI95.1"]> 1 & rr[x, "CI95.2"]> 1, 2, 0))  })
  
  ci <- sapply(1:length(rr$cause), function(x){paste(round(rr[x,"CI95.1"], 3), 
                                                     round(rr[x, "CI95.2"], 3),
                                                     sep= "-") })
  cbind(t(alpha), t(rr$RR),  t(ci), t(rr$p.value))
  
} # return 3 number corresponding to significance of ACM rr, OCM rr, PCSM rr


## function that from RR conf interval tells if is significant or not (1 = risk reduction, 0 = no difference, 2 = higher risk)
alpha_sign_v3 <- function (rr){
  # 1 for significant  reduction result, 0 for non significant result or higher risk for invited
  alpha <- sapply(1:length(rr$cause), function(x){ifelse(rr[x,"CI95.1"]< 1 & rr[x, "CI95.2"]< 1, 1, 
                                                         ifelse(rr[x,"CI95.1"]> 1 & rr[x, "CI95.2"]> 1, 2, 0))  })
  
  ci <- sapply(1:length(rr$cause), function(x){paste(round(rr[x,"CI95.1"], 3), 
                                                     round(rr[x, "CI95.2"], 3),
                                                     sep= "-") })
  cbind(t(alpha), (t(rr$R.p0)), 
                  (t(rr$R.p1)), 
                  (t(rr$RR)), 
                  (t(rr$SE_log.RR)), 
                  (t(rr$Z)),  t(ci), 
                  (t(rr$p.value)))
  
} # return 3 number corresponding to significance of ACM rr, OCM rr, PCSM rr



RR.result <- function(i, sfit, time) {
  
  if(i == 1){
    pcif <- 1 - summary(sfit, time = time)$pstate[, i]
  } else {
    pcif <- summary(sfit, time = time)$pstate[, i]
  }
  
  secif <- summary(sfit, time = time)$std.err[, i]
  names(secif) <- names(pcif) <- c("p0", "p1")
  V <- diag(secif^2)
  
  # 95% normal-based delta-method CI for logRR + exp
  G <- matrix(
    c(
      -1/pcif["p0"], # deriv wrt p0
      1/pcif["p1"] # deriv wrt p1
    ),
    ncol = 2) # Gradient
  ase <- sqrt(G %*% V %*% t(G))  ## SE(log RR)
  
  X <- c(
    cause = i, 
    time =  time, 
    R = pcif["p0"], 
    R = pcif["p1"], 
    RR = as.numeric(pcif["p1"]/pcif["p0"]),
    SE_log.RR = ase,
    Z = log(pcif["p1"]/pcif["p0"])/ase, 
    p.value = 2*pnorm(-abs(log(pcif["p1"]/pcif["p0"])/ase), lower.tail = TRUE),
    CI95. = (exp(log(pcif["p1"]/pcif["p0"]) + c(-1, 1)*qnorm(.975)*c(ase))), 
    CI98. = (exp(log(pcif["p1"]/pcif["p0"]) + c(-1, 1)*qnorm(.99)*c(ase)))    ) 
  return(X)
  
}


surv.analysis.RR <- function(dataset) {
  sfit <- survfit(Surv(time.predicted, event.predicted) ~ group_rand, data = dataset)
  time.RR <- max(sfit[["time"]])
  RR <- rbind(
    RR.result(1, sfit = sfit, time = time.RR),  ## 1 - survival
    RR.result(2, sfit = sfit, time = time.RR),  ## other cause
    RR.result(3, sfit = sfit, time = time.RR)   ## prostate cancer
  ) %>% as.data.frame()
  RR$cause <- c("all.cause", "other.cause", "pca")
  
  output <- 
    c(table(dataset$event.predicted),
      alpha_sign_v2(RR) )
  # alpha sign -> 0 non significant, 1 reduction, 2 increase
  
  return(list (Surv = sfit , RR.an = output, risks = RR))
}

surv.analysis.RR3 <- function(dataset) {
  sfit <- survfit(Surv(time.predicted, event.predicted) ~ group_rand, data = dataset)
  time.RR <- max(sfit[["time"]])
  RR <- rbind(
    RR.result(1, sfit = sfit, time = time.RR),  ## 1 - survival
    RR.result(2, sfit = sfit, time = time.RR),  ## other cause
    RR.result(3, sfit = sfit, time = time.RR)   ## prostate cancer
  ) %>% as.data.frame()
  RR$cause <- c("all.cause", "other.cause", "pca")
  
  output <- 
    c(table(dataset$event.predicted),
      alpha_sign_v3(RR) )
  names(output) <- c("Alive", "Other death", "PC death",
                     ## from alpha_sign function   # cbind(t(alpha),t(rr$R0), t(rr$R1),  t(rr$RR), t(rr$SE_log.RR), t(rr$Z),  t(ci), t(rr$p.value))
                     "a_ACM", "a_OCM", "a_PCSM",
                     "R0_ACM", "R0_OCM", "R0_PCSM",
                     "R1_ACM", "R1_OCM", "R1_PCSM", 
                     "RR_ACM", "RR_OCM", "RR_PCSM", 
                     "ase_ACM", "ase_OCM", "ase_PCSM",
                     "Z_ACM", "Z_OCM", "Z_PCSM",
                     "CI_ACM", "CI_OCM", "CI_PCSM", 
                     "p.value_ACM", "p.value_OCM", "p.value_PCSM")    
  # alpha sign -> 0 non significant, 1 reduction, 2 increase
  
  return(list (Surv = sfit , RR.an = output, risks = RR))
}


