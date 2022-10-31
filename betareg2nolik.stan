data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> Narea;
  int<lower=0> T;
  int iiobs[Narea,T];
  real<lower=0,upper=1> y[N];
  matrix[N, K] X;
  matrix[Narea,Narea] cormat;
}

transformed data {
  matrix[Narea,Narea] L;
  real<lower=0,upper=1> yvec[T,Narea];
  for (i in 1:Narea) {
    for (t in 1:T) {
       yvec[t][i] = y[(t-1)*Narea+i];
    }
  }
  L = cholesky_decompose(cormat);
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
     omega ~ multi_normal_cholesky(rep_vector(0,Narea),L);
}

generated quantities {
  vector[N] log_lik;  
  {
  vector[N] xb;
  vector[N] O;
  vector[N] iiobs2;
  xb = X * betas;
  for (i in 1:Narea) {
    for (t in 1:T) {
      O[(t-1)*Narea+i] = sigma * omega[i];
      iiobs2[(t-1)*Narea+i] = iiobs[i,t];
  }
  }
  for (i in 1:N) if (iiobs2[i] > 0) log_lik[i] = beta_proportion_lpdf(y[i] | inv_logit(xb[i] + O[i]), kappa);
  }
}