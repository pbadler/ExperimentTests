data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;                // observation vector
  vector[N] X;                // size vector
  int<lower=0> Wcovs;
  matrix[N, Wcovs] W;
}
parameters{
  // for training data model  
  real<lower=0> sigma;
  real<lower=0>b1;
  real<lower=0>b2;
  vector[Wcovs] w; 
}
transformed parameters{
  vector[N] mu;
  vector[N] a; 
  
  a <- W*w;
  
  for( n in 1:N){
    mu[n] <- log( a[n]*exp(-b1*b2^X[n]));
  }  
}
model{
   // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  a ~ normal(0,10);
  b1 ~ normal(0,10);
  b2 ~ normal(0,10);
  w ~ normal(0, 10); 
  
  // Likelihood
  Y ~ normal(mu, sigma);
  
}
generated quantities {
  vector[N] log_lik;          // vector for computing log pointwise predictive density  
  vector[N] y_hat;
  
  for(n in 1:N){
      y_hat[n] <- normal_rng(mu[n], sigma);
      log_lik[n] <- normal_log(Y[n], mu[n], sigma);
  }
}

