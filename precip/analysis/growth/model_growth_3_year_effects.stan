data{
  // training datalist, historical observations 
  int<lower=0> N2;             // observations
  vector[N2] Y2;  // observation vector
  vector[N2] X2;                // size vector
  int<lower=0> G;             // groups
  matrix[N2, G] gm2;            // group dummy variable matrix 
  int<lower=0> nyrs2;          // years
  int<lower=0> yid2[N2];        // year id
  int<lower=0> Wcovs;         // number of crowding effects 
  matrix[N2,Wcovs] W2;          // crowding matrix
  real<lower=0> tau_beta;    // prior sd for climate effects 
}
parameters{
  // for training data model  
  vector[G] bg;                     // varying group effects with first group as intercept 
  vector[nyrs2] a_raw;
  real b1_mu;
  vector[nyrs2] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[Wcovs] w;
  real<lower=0> tau; 
  real<lower=0> tauSize; 
}
transformed parameters{
// for training data model  
  vector[nyrs2] a;
  vector[nyrs2] b1;
  vector[N2] gint; 
  real mu[N2];
  vector[N2] crowdEff;
  vector[N2] climEff; 
  vector[N2] sigma;

  // for training data model -----------------------------------
  gint <- gm2*bg;
  crowdEff <- W2*w;

  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N2){
    mu[n] <- gint[n] + a[yid2[n]] + b1[yid2[n]]*X2[n] + crowdEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
  }
}
model{
   // for training data model 
  // Priors
  //sigma ~ cauchy(0, 5);
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  tau ~ normal( 0, 1); 
  tauSize ~ normal(0,1);
  w ~ normal(0,tau_beta);

  // Likelihood
  Y2 ~ normal(mu, sigma);
}
generated quantities {
  // hold out predictions 
  vector[N2] log_lik;          // vector for computing log pointwise predictive density  

  for(n in 1:N2){
      log_lik[n] <- normal_log(Y2[n], mu[n], sigma[n]);
  }
}
