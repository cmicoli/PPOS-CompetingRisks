## Functions for PPOS analysis in the ISPY-COVID example

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


# Function to fit weibull models (PHASE 1 - modelling), stratified over treatment arm (no PH assumption)
fit.noPH <- function(dataset, 
                     brms_function1, 
                     brms_function2
){  ### models fitted within treatment arm
  fit1 <- brm(brms_function1, data = dataset,     
              chains = 4,
              iter = 3500,
              warmup = 500,
              seed = 12345,
              cores = 4,
              stanvars = stanvars_weibullPH ,
              prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                        prior(normal(0, sqrt(.5)), class = b), 
                        prior(exponential(1), class = gamma)),
              backend = "cmdstanr")
  
  # model for DEATH before recovery
  fit2 <- brm(brms_function2, data = dataset,
              chains = 4,
              iter = 3500,
              warmup = 500,
              seed = 12345,
              cores = 4,
              stanvars = stanvars_weibullPH,
              prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                        prior(normal(0, sqrt(.5)), class = b), 
                        prior(exponential(1), class = gamma)),
              backend = "cmdstanr")
  
  return(list(Fit1 = fit1, Fit2= fit2 ) )
}


# Function to fit weibull models (PHASE 1 - modelling) with prior on treatment effect (PH assumption on treatment effect)
fit.PH.HR <- function(dataset, 
                      brms_function1, 
                      brms_function2, 
                      HR.mean , HR.sd){  ### prior for HR in mortality is ->  normal(log(1), sqrt(0.5)) 
  stanvars <- stanvar(HR.mean, name='mean') + stanvar(HR.sd, name='sd')
  fit1 <- brm(brms_function1, data = dataset,     
              chains = 4,
              iter = 3500,
              warmup = 500,
              seed = 12345,
              cores = 4,
              stanvars = stanvars_weibullPH + stanvars,
              prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                        prior(normal(mean, sd), class = "b", coef = "TREATMENTNEW_TR"),
                        prior(normal(0, sqrt(.5)), class = b), 
                        prior(exponential(1), class = gamma)),
              backend = "cmdstanr")
  
  # model for DEATH before recovery
  fit2 <- brm(brms_function2, data = dataset,
              chains = 4,
              iter = 3500,
              warmup = 500,
              seed = 12345,
              cores = 4,
              stanvars = stanvars_weibullPH,
              prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                        prior(normal(log(1), sqrt(0.5)), class = "b", coef = "TREATMENTNEW_TR"),
                        prior(normal(0, sqrt(.5)), class = b), 
                        prior(exponential(1), class = gamma)),
              backend = "cmdstanr")
  
  return(list(Fit1 = fit1, Fit2= fit2, model = c(mean = HR.mean, sd = HR.sd)) )
}


# Function to predict data (PHASE 2 - prediction) - when fitting models with no PH assumption
new_predicted_dataset_noPH <- function (data_to_predict, data_observed, n.sim, 
                                        model1, model2 , cores){
  
  r.beta.HR <- random.draw.beta.weibull(model1, model2, n.sim)
  
  predicted <- mclapply( 1:n.sim, 
                         function(x) {
                           predicted <- rand_crisk_w(
                             dim(data_to_predict[[x]])[1],
                             as.numeric(r.beta.HR$fit.1[x,]),
                             as.numeric(r.beta.HR$fit.2[x,]),
                             model.matrix(~ SCRNCOVID_bin, data = data_to_predict[[x]]), ### checkkk!!!!
                             data_to_predict[[x]]$rec_FU ### dat_inv0$maxFU18??? or time of real censoring
                           )
                           
                           predicted <- 
                             cbind(data_to_predict[[x]] , predicted)
                           
                           predicted_cens <- predicted
                           predicted_cens$time.predicted <- ifelse(predicted$time <= 60, predicted$time, 60)
                           predicted_cens$event.predicted <- factor(ifelse(predicted$time <= 60, predicted$event, "0"))
                           
                           predicted_dat_tot <- 
                             data_observed %>% mutate(event = NA, 
                                                      time  = NA, 
                                                      event.predicted = EVENT, 
                                                      time.predicted = rec_FU) %>%
                             #rbind(predicted)
                             rbind(predicted_cens)
                           
                           predicted_dat_tot
                           
                         },
                         mc.cores = cores)
  
  return(predicted)
}


# Function to predict data (PHASE 2 - prediction) - when fitting models with PH assumption
predA <- function (data_to_predict, data_observed, n.sim, 
                   model1, model2, cores ){  ### with PH for treatment
  
  r.beta.HR <- random.draw.beta.weibull(model1, model2, n.sim)
  #print(r.beta.HR)
  predicted <- mclapply( 1:n.sim, 
                         function(x) {
                           predicted <- rand_crisk_w(
                             dim(data_to_predict[[x]])[1],
                             as.numeric(r.beta.HR$fit.1[x,]),
                             as.numeric(r.beta.HR$fit.2[x,]),
                             model.matrix(~ TREATMENT + SCRNCOVID_bin, data = data_to_predict[[x]]), ### checkkk!!!!
                             data_to_predict[[x]]$rec_FU ### dat_inv0$maxFU18??? or time of real censoring
                           )
                           
                           predicted <- 
                             cbind(data_to_predict[[x]] , predicted)
                           
                           predicted_cens <- predicted
                           predicted_cens$time.predicted <- ifelse(predicted$time <= 60, predicted$time, 60)
                           predicted_cens$event.predicted <- factor(ifelse(predicted$time <= 60, predicted$event, "0"))
                           
                           predicted_dat_tot <- 
                             data_observed %>% mutate(event = NA, 
                                                      time  = NA, 
                                                      event.predicted = EVENT, 
                                                      time.predicted = rec_FU) %>%
                             #rbind(predicted)
                             rbind(predicted_cens)
                           
                           predicted_dat_tot
                           
                         },
                         mc.cores = cores)
}


# compiling model for analysis, than will then be update. Main analysis in the trial to evaluate success
compile.analysis.model <- function(dataset) {
  dataset$event.all <- ifelse(dataset$event.predicted %in% c(1,2), 1, 0)
  dataset$event1 <- dataset$event.all * (dataset$event.predicted == 1)  ## "Recovery"
  dataset$event2 <- dataset$event.all * (dataset$event.predicted == 2)  ## "Death w/o recovery"
  dataset$censored1 <- 1 - dataset$event1
  dataset$censored2 <- 1 - dataset$event2
  
  formula <- bf(time.predicted | cens(censored1) ~   0 + Intercept + TREATMENT + SCRNCOVID_bin,
     family = weibullPH)
  fit.null <- brm(formula,
                  data = dataset,     
                  chains = 0,
                  cores = 1,
                  stanvars = stanvars_weibullPH,
                  prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                            prior(normal(0, sqrt(.5)), class = b), 
                            prior(exponential(1), class = gamma)),
                  backend = "cmdstanr")
  fit.null
}


# Main analysis in the trial to evaluate success, done over the "dataset", and using the compired model
main.analysis.update <- function(dataset, fit.null){ ## event considered is conteined in EVENT.PREDICTED
  
  ## starting from main dataset and formula to analyize data 
  ## setting up the dataset with the proper variable
  dataset$event.all <- ifelse(dataset$event.predicted %in% c(1,2), 1, 0)
  dataset$event1 <- dataset$event.all * (dataset$event.predicted == 1)  ## "Recovery"
  dataset$event2 <- dataset$event.all * (dataset$event.predicted == 2)  ## "Death w/o recovery"
  dataset$censored1 <- 1 - dataset$event1
  dataset$censored2 <- 1 - dataset$event2
  ## output is median HR, and 95% qCrI for treatment effect, 
  ## and probability of recovery graduating : p(HR >1 | Dati) = # draws with HR >1 / 4000
  fit.main <- update(fit.null,
                     recompile = FALSE,
                     newdata = dataset,     
                     chains = 4,
                     iter = 1500,
                     warmup = 500,
                     seed = 12345,
                     cores = 1,
                     stanvars = stanvars_weibullPH,
                     prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                               prior(normal(0, sqrt(.5)), class = b), 
                               prior(exponential(1), class = gamma)),
                     backend = "cmdstanr") ## rstan
  #fit.main
  pos.draws <- fit.main %>% 
    tidybayes::spread_draws(b_Intercept, b_TREATMENTNEW_TR , b_SCRNCOVID_bin6D7   , gamma)
  # print(dim(pos.draws))
  # print(colSums(is.na(pos.draws)))
  HR.CrI <- quantile(exp(pos.draws$b_TREATMENTNEW_TR),
                     probs =c(0.5, 0.025, 0.975),  na.rm = FALSE) %>% 
    as.data.frame() %>% t() 
  HR.CrI.pr <- HR.CrI %>% 
    cbind(count.HRvs1  = sum(exp(pos.draws$b_TREATMENTNEW_TR)>1)) %>% 
    as.data.frame() %>%
    mutate(pr = count.HRvs1 /4000, ## probability of graduating in this dataset
           success = ifelse(pr > 0.975, 1, 0 )) ## definition of success with this dataset
  colnames(HR.CrI.pr) <- c( "HR.median",  "HR.025lower" ,  "HR.975upper"   ,   "count.HRvs1",  "pr", "success")
  return(list(post.draws = pos.draws, analysis = HR.CrI.pr))
}




main.analysis <- function(dataset, formula, n.cores){ ## event considred is conteined in EVENT.PREDICTED
  
  ## starting from main dataset and formula to analyize data 
  ## setting up the dataset with the proper variable
  dataset$event.all <- ifelse(dataset$event.predicted %in% c(1,2), 1, 0)
  dataset$event1 <- dataset$event.all * (dataset$event.predicted == 1)  ## "Recovery"
  dataset$event2 <- dataset$event.all * (dataset$event.predicted == 2)  ## "Death w/o recovery"
  dataset$censored1 <- 1 - dataset$event1
  dataset$censored2 <- 1 - dataset$event2
  ## output is median HR, and 95% qCrI for treatment effect, 
  ## and probability of recovery graduating : p(HR >1 | Dati) = # draws with HR >1 / 4000
  fit.main <- brm(formula,
                  data = dataset,     
                  chains = 4,
                  iter = 1500,
                  warmup = 500,
                  seed = 12345,
                  cores = n.cores,
                  stanvars = stanvars_weibullPH,
                  prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                            prior(normal(0, sqrt(.5)), class = b), 
                            prior(exponential(1), class = gamma)),
                  backend = "cmdstanr") ## rstan
  fit.main
  pos.draws <- fit.main %>% 
    tidybayes::spread_draws(b_Intercept, b_TREATMENTNEW_TR , b_SCRNCOVID_bin6D7   , gamma)
  # print(dim(pos.draws))
  # print(colSums(is.na(pos.draws)))
  HR.CrI <- quantile(exp(pos.draws$b_TREATMENTNEW_TR),
                     probs =c(0.5, 0.025, 0.975),  na.rm = FALSE) %>% 
    as.data.frame() %>% t() 
  HR.CrI.pr <- HR.CrI %>% 
    cbind(count.HRvs1  = sum(exp(pos.draws$b_TREATMENTNEW_TR)>1)) %>% 
    as.data.frame() %>%
    mutate(pr = count.HRvs1 /4000, ## probability of graduating in this dataset
           success = ifelse(pr > 0.975, 1, 0 )) ## definition of success with this dataset
  return(list(post.draws = pos.draws, analysis = HR.CrI.pr))
}


# Main analysis for success, to be runned on interim data
interim.main.analysis <- function(dataset, formula1, formula2, n.cores){
  ## event considred is conteined in EVENT
  
  ## starting from main dataset and formula to analyize data 
  ## setting up the dataset with the proper variable
  dataset$event.all <- ifelse(dataset$EVENT %in% c(1,2), 1, 0)
  dataset$event1 <- dataset$event.all * (dataset$EVENT == 1)  ## "Recovery"
  dataset$event2 <- dataset$event.all * (dataset$EVENT == 2)  ## "Death w/o recovery"
  dataset$censored1 <- 1 - dataset$event1
  dataset$censored2 <- 1 - dataset$event2
  ## output is median HR, and 95% qCrI for treatment effect, 
  ## and probability of recovery graduating : p(HR >1 | Dati) = # draws with HR >1 / 4000
  fit.main.event1 <- brm(formula1,
                  data = dataset,     
                  chains = 4,
                  iter = 1500,
                  warmup = 500,
                  seed = 12345,
                  cores = n.cores,
                  stanvars = stanvars_weibullPH,
                  prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                            prior(normal(0, sqrt(.5)), class = b), 
                            prior(exponential(1), class = gamma)),
                  backend = "cmdstanr") ## rstan
  fit.main.event2 <- brm(formula2,
                  data = dataset,     
                  chains = 4,
                  iter = 1500,
                  warmup = 500,
                  seed = 12345,
                  cores = n.cores,
                  stanvars = stanvars_weibullPH,
                  prior = c(prior(normal(0, 20), class = "b", coef = "Intercept"),
                            prior(normal(0, sqrt(.5)), class = b), 
                            prior(exponential(1), class = gamma)),
                  backend = "cmdstanr") ## rstan
  fit.main.event1
  pos.draws <- fit.main.event1 %>% 
    tidybayes::spread_draws(b_Intercept, b_TREATMENTNEW_TR , b_SCRNCOVID_bin6D7   , gamma)
  # print(dim(pos.draws))
  # print(colSums(is.na(pos.draws)))
  HR.CrI <- quantile(exp(pos.draws$b_TREATMENTNEW_TR),
                     probs =c(0.5, 0.025, 0.975),  na.rm = FALSE) %>% 
    as.data.frame() %>% t() 
  HR.CrI.pr <- HR.CrI %>% 
    cbind(count.HRvs1  = sum(exp(pos.draws$b_TREATMENTNEW_TR)>1)) %>% 
    as.data.frame() %>%
    mutate(pr = count.HRvs1 /4000, ## probability of graduating in this dataset
           success = ifelse(pr > 0.975, 1, 0 )) ## definition of success with this dataset
  return(list(fit1 = fit.main.event1, fit2 = fit.main.event2, 
              post.draws = pos.draws, analysis = HR.CrI.pr))
}








