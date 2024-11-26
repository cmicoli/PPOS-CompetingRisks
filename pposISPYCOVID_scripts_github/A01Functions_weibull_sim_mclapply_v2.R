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

