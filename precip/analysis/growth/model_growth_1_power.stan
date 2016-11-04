data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  vector[N] X;                // size vector
}
parameters{
  // for training data model  
  real<lower=0> sigma; 
  real<lower=0>K;
  real b;
}
transformed parameters{
  vector[N] mu;
  
  for( n in 1:N){
    mu[n] <- (X[n]*K)/(X[n]+(K-X[n])*exp(b));
  }  
}
model{
   // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  K ~ normal(0,10);
  b ~ normal(0,10);

  // Likelihood
  Y ~ normal(mu, sigma);
  
}
generated quantities {
  vector[N] log_lik;          // vector for computing log pointwise predictive density  

  for(n in 1:N){
      log_lik[n] <- normal_log(Y[n], mu[n], sigma);
  }
}

