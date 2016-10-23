data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;                // observation vector
  vector[N] X;                // size vector
  int<lower=0> G;             // groups
  matrix[N, G] gm;            // group dummy variable matrix 
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  vector[N] W;                // intraspecific crowding
}
parameters{
  // for training data model  
  vector[G] bg;                     // varying group effects with first group as intercept 
  real b1_mu;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[nyrs] a_raw;
  vector[nyrs] b1_raw;
  real tauSize;
  real tau;
  real w;
}
transformed parameters{
  vector[nyrs] a;
  vector[nyrs] b1;
  
  // for training data model  
  real<lower=0> sigma[N];
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;

  // for training data model -----------------------------------
  gint     <- gm*bg;
  crowdEff <- W*w;
  
  // reparamaterize the hierarchical parameters  
  a <- 0 + sig_a*a_raw;
  b1 <- b1_mu + sig_b1*b1_raw;

  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
  }

}
model{
  // for training data model 
  // Priors
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  
  //a ~ normal(0, sig_a); 
  //b1 ~ normal(b1_mu, sig_b1);
  tau ~ normal(0,10);
  tauSize ~ normal(0,10);
  w ~ normal(0,10);

  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {
  // hold out predictions 
  vector[N] log_lik;                      // vector for computing log pointwise predictive density  

  for(n in 1:N){
      log_lik[n] <- normal_log(Y[n], mu[n], sigma[n]);
  }
}
