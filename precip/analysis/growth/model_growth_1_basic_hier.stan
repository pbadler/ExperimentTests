data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  vector[N] X;                // size vector
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
}
parameters{
  // for training data model  
  real<lower=0> sigma; 
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real alpha; 
}
transformed parameters{
// for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  real mu[N];

  // for training data model -----------------------------------
  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- alpha + a[yid[n]] + b1[yid[n]]*X[n];
  }
  
}
model{
   // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  alpha ~ normal(0,10);

  // Likelihood
  Y ~ normal(mu, sigma);
  
}
generated quantities {
  // hold out predictions 
  vector[N] log_lik;          // vector for computing log pointwise predictive density  

  for(n in 1:N){
      log_lik[n] <- normal_log(Y[n], mu[n], sigma);
  }
}

