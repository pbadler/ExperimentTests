data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;                // observation vector
  vector[N] X;                // size vector
}
parameters{
  // for training data model  
  real<lower=0> sigma;
  real<lower=0>K;
  real<lower=0>r;
}
transformed parameters{
  vector[N] mu;
  vector[N] y_log; 
  
  for( n in 1:N){
    mu[n] <- log(K*(X[n]/K)^exp(-r));
    y_log[n] <- log(Y[n]);
  }  
}
model{
  // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  K ~ normal(0,10);
  r ~ normal(0,10);
  
  // Likelihood
  y_log ~ normal(mu, sigma);
  
}
generated quantities {
  vector[N] log_lik;          // vector for computing log pointwise predictive density  

  for(n in 1:N){
      log_lik[n] <- normal_log(y_log[n], mu[n], sigma);
  }
}

