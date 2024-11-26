cdf_all_w <- function(t, beta1, beta2, model_matrix) {
  stopifnot(length(beta1) == length(beta2))
  #stopifnot(length(t) == 1)
  
  shape_coeff_position <- length(beta1)
  
  1 - exp(
    -(
      flexsurv::HweibullPH(t, 
                           scale = exp(model_matrix %*% beta1[-shape_coeff_position]), 
                           shape = beta1[shape_coeff_position]) +  
        flexsurv::HweibullPH(t, 
                             scale = exp(model_matrix %*% beta2[-shape_coeff_position]), 
                             shape = beta2[shape_coeff_position]) 
    )
  )
}


cdf_all_trunc_w <- function(x, beta1, beta2, model_matrix, trunc) {
  cdf_at_truncation <- ifelse(trunc <= 0, 0, cdf_all_w(t = trunc, beta1, beta2, model_matrix))
  
  ifelse(x < trunc, 0, (cdf_all_w(t = x, beta1, beta2, model_matrix) - cdf_at_truncation) / (1 - cdf_at_truncation))
}



quantilef_all_w <- function(q, beta1, beta2, model_matrix, trunc) {
  sapply(
    seq_along(q), function(x) {
      uniroot(function(z) {
        cdf_all_trunc_w(z, beta1 = beta1, beta2 = beta2, model_matrix = model_matrix[x, ], trunc[x]) - q[x]  
      }, lower = 0, upper  = 400, extendInt = "upX")$root
    }
  )
}



rand_crisk_w <- function(n, beta1, beta2, model_matrix, trunc) {
  # uniform random values
  u <- runif(n, min = 0, max = 1)
  
  # step 2 of Beyersmann (from u to random event times - we don't know what competing event happened yet)
  event_times <- quantilef_all_w(u, beta1 = beta1, beta2 = beta2, model_matrix = model_matrix, trunc)
  
  # step 3 of Beyersmann (what competing event happened at each event time?)
  shape_coeff_position <- length(beta1)
  h_1 <- flexsurv::hweibullPH(event_times, 
                              scale = exp(model_matrix %*% beta1[-shape_coeff_position]), 
                              shape = beta1[shape_coeff_position]) 
  h_2 <- flexsurv::hweibullPH(event_times, 
                              scale = exp(model_matrix %*% beta2[-shape_coeff_position]), 
                              shape = beta2[shape_coeff_position]) 
  
  causes <- rbinom(n, 1, prob = h_2/(h_1+h_2)) + 1  
  
  data.frame(
    time = as.numeric(event_times),
    event = as.numeric(causes)  
  )
}


random.draw.beta.weibull <- function(fit1, fit2, n.sim) {
  length.fit1 <- dim(as.data.frame(fit1))[2]
  length.fit2 <- dim(as.data.frame(fit2))[2]
  fit <- list(
    fit.1 = (as.data.frame(fit1)[, c(-length.fit1, -length.fit1+1)]),
    fit.2 = (as.data.frame(fit2)[, c(-length.fit2, -length.fit2+1)])
  )
  
  #r.beta <- lapply(fit, function(x) {x[sample(nrow(x), 1),]})
  r.beta <- lapply(fit, function(x) {x[do.call("sample", list(nrow(x),n.sim)),]})  
  #r.beta <- lapply(r.beta, function(x) {-as.numeric(x)})
  
  r.beta
  
}





