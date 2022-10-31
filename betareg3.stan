data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> Narea;
  int<lower=0> T;
  int iiobs[Narea,T];
  real<lower=0,upper=1> y[N];
  matrix[N, K] X;
  matrix[Narea,Narea] cormat;
  int<lower=0> Nobs;
  int<lower=0> iobsvec[N];
}

transformed data {
  matrix[Narea,Narea] L;
  vector[Narea] rep0;
  real<lower=0,upper=1> yvec[T,Narea];
  int iiobs2[N];
  for (i in 1:Narea) {
    for (t in 1:T) {
       yvec[t][i] = y[(t-1)*Narea+i];
       iiobs2[(t-1)*Narea+i] = iiobs[i,t];
    }
  }
  L = cholesky_decompose(cormat);
  rep0 = rep_vector(0,Narea)
}

parameters {
  vector[K] betas;
  real<lower=0.001> kappa;
  vector[Narea] omega;
  real<lower=0> sigma;
}

model {
     for (t in 1:T) {
        for (i in 1:Narea) {
          if (iiobs[i,t] > 0) yvec[t][i] ~ beta_proportion(inv_logit(X[(t-1)*Narea+i,] * betas + sigma * omega),kappa);	 
        }
     }
     kappa ~ gamma(2,1);
     sigma ~ gamma(2,1);    
     omega ~ multi_normal_cholesky(rep0,L);
}

generated quantities {
  vector[Nobs] log_lik;  
  {
  vector[N] xb;
  vector[N] O;
  xb = X * betas;
  for (i in 1:Narea) {
    for (t in 1:T) {
      O[(t-1)*Narea+i] = sigma * omega[i];
  }
  }
  for (i in 1:N) if (iiobs2[i] > 0) log_lik[iobsvec[i]] = beta_proportion_lpdf(y[i] | inv_logit(xb[i] + O[i]), kappa);   
  }
}