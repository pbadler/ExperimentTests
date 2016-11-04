data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  vector[N] X;                // size vector
}
parameters{
  //real<lower=0> sigma; 
  real b0;
  real b1;
  real tau;
  real tauSize;
}
transformed parameters{
  vector[N] mu;
  vector[N] X_log; 
  vector[N] Y_log;
  vector[N] sigma;
  
  for( n in 1:N){
    X_log[n] <- log(X[n]);
    Y_log[n] <- log(Y[n]); 
  } 
  
  mu <- b0 + b1*X_log;
  
  for( n in 1:N)
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  

}
model{
   // for training data model 
  // Priors
  //sigma ~ cauchy(0, 5);
  b0 ~ normal(0,10);
  b1 ~ normal(0,10);
  tau ~ normal(0,10);
  tauSize ~ normal(0,10);
  
  // Likelihood

  Y_log ~ normal(mu, sigma);
  
}
generated quantities {
  vector[N] log_lik;          // vector for computing log pointwise predictive density  

  for(n in 1:N){
      log_lik[n] <- normal_log(Y_log[n], mu[n], sigma);
  }
}

