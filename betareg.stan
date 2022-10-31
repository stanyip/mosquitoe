data {
  int<lower=0> N;
  int<lower=0> K;
  real<lower=0,upper=1> y[N];
  matrix[N, K] X;
}

parameters {
  vector[K] betas;
  real<lower=0.001> kappa;
}

model {
     y ~ beta_proportion(inv_logit(X * betas),kappa);
     kappa ~ gamma(2,1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] xb;
  xb = X * betas;
  for (i in 1:N)
    log_lik[i] = beta_proportion_lpdf(y[i] | inv_logit(xb[i]), kappa);
}