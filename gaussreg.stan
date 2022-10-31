data {
  int<lower=0> N;
  int<lower=0> K;
  real<lower=0,upper=1> y[N];
  matrix[N, K] X;
}

parameters {
  vector[K] betas;
  real<lower=0.001> sigma;
}

model {
     y ~ normal(X * betas,sigma);
     sigma ~ gamma(2,1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] xb;
  xb = X * betas;
  for (i in 1:N)
    log_lik[i] = normal_lpdf(y[i] | xb[i], sigma);
}