data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> Narea;
  int<lower=0> T;
  real<lower=0,upper=1> y[N];
  matrix[N, K] X;
  matrix[Narea,Narea] cormat;
  int<lower=1,upper=Narea> areano[N];
  int<lower=1,upper=12> month[N];
}

transformed data {
  matrix[Narea,Narea] L;
  vector[Narea] rep0;
  vector[N] summer;
  L = cholesky_decompose(cormat);
  rep0 = rep_vector(0,Narea);
  for (i in 1:N) {
      if ((month[i] >=6) && (month[i] <=8)) summer[i] = 1; else summer[i] = 0;
  }
}

parameters {
  vector[K] betas;
  real<lower=0.001> kappa;
  vector[Narea] omega;
  vector[N] epsilon;
  real<lower=0> sigma;
  real<lower=0> sigmao;
}

transformed parameters {
vector[N] mu;
     for (i in 1:N) {
          mu[i] = inv_logit(X[i,] * betas + summer[i] * sigmao * omega[areano[i]] + sigma * epsilon[i]);
     }
}

model {
     y ~ beta_proportion(mu,kappa);	 
     kappa ~ gamma(2,1);
     sigmao ~ gamma(2,1);    
     sigmao ~ gamma(2,1);    
     omega ~ multi_normal_cholesky(rep0,L);
     epsilon ~ normal(0,1);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = beta_proportion_lpdf(y[i] | mu[i], kappa);
}
