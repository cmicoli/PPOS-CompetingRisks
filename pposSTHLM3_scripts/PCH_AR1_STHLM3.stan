// PCH AR1 STHLM3
// generated with brms 2.20.12
// edited by Chiara & Andrea 20241009
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  array[N] int<lower=-1,upper=2> cens;  // indicates censoring
  array[N] real lb;  // lower truncation bounds;
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
parameters {
  vector[K] b;  // regression coefficients
  real<lower=0> tau;
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b[1] | 10, 20);
  for (i in 2:(K-1)) {
    lprior += normal_lpdf(b[i] - b[i-1] | 0, tau);
  }
  lprior += normal_lpdf(b[K] | 0, sqrt(.5));
  lprior += exponential_lpdf(tau | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += X * b;
    mu = exp(mu);
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        target += exponential_lpdf(Y[n] | inv(mu[n])) -
        exponential_lccdf(lb[n] | inv(mu[n]));
      } else if (cens[n] == 1) {
        target += exponential_lccdf(Y[n] | inv(mu[n])) -
        exponential_lccdf(lb[n] | inv(mu[n]));
      } else if (cens[n] == -1) {
        target += exponential_lcdf(Y[n] | inv(mu[n])) -
        exponential_lccdf(lb[n] | inv(mu[n]));
      }
    }
  }
    // priors including constants
  target += lprior;
}
generated quantities {
}
