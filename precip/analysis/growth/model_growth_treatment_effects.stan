data{
  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  vector[Nhold] Yhold;
  int<lower=0> nyrs;                // training data years 
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      // year out id
  int<lower=0> G;                   // Groups 
  matrix[Nhold, G] gmhold;          // group dummy variable matrix
  int<lower=0> nThold;              // number of treatments 
  matrix[Nhold, nThold] tmhold;     // treatment dummy matrix  
  int<lower=0> Wcovs;               // number of crowding effects 
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
}
parameters{
  // for training data model  
  vector[nyrshold] a_raw;
  real b1_mu;
  vector[nyrshold] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sigma;
  vector[Wcovs] w;
  vector[nThold] bt;
  vector[G] bg;
}
transformed parameters{
  // for training data model  
  vector[nyrshold] a;
  vector[nyrshold] b1;
  vector[Nhold] treatEff;
  vector[Nhold] gint;
  real mu[Nhold];
  vector[Nhold] crowdEff;
  vector[Nhold] year_effect; 
  
  // for training data model -----------------------------------
  crowdEff <- Whold*w;
  treatEff <- tmhold*bt;

  gint <- gmhold*bg; 
  
  a  <- 0 + sig_a*a_raw; 
  b1 <- 0 + sig_b1*b1_raw; 
    
  for(n in 1:Nhold){

    year_effect[n] <- a[yidhold[n] - nyrs] + b1[yidhold[n] - nyrs]*Xhold[n];
    
    mu[n] <- gint[n] + treatEff[n] + year_effect[n] + b1_mu*Xhold[n] + crowdEff[n];
  }
  
}
model{
  // Priors
  bt ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);

  // Likelihood
  Yhold ~ normal(mu, sigma);
}
generated quantities {
  
  vector[Nhold] log_lik2;                     // for heldout data 

  # fitted data log_lik 
  for(n in 1:Nhold){
      log_lik2[n] <- normal_log(Yhold[n], mu[n], sigma);
  }
}
