data{
  int<lower=0> N;             // observations
  int<lower=0,upper=1> Y[N];  // observation vector
  vector[N] X;                // size vector
  int<lower=0> G;             // groups
  matrix[N, G] gm;            // group dummy variable matrix 
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  int<lower=0> Wcovs;         // number of crowding effects 
  matrix[N,Wcovs] W;          // crowding matrix
  int<lower=0> spp;           // focal species 
  
   // for out of sample predictions
  int<lower=0> Nhold;
  int<lower=0,upper=1> Yhold[Nhold];
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  int<lower=0> gidhold[Nhold];      // group id holdout
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data

  
}
parameters{
  vector[G] bg;                     // varying group effects with first group as intercept 
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real w;
}
transformed parameters{
  vector[nyrs] a;
  vector[nyrs] b1;
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] W_intra; 

  // group, crowding and climate effects 
  gint <- gm*bg;
  W_intra <- W[, spp];
  crowdEff <- W_intra*w;
  
  // reparamaterize the hierarchical parameters  
  a <- 0 + sig_a*a_raw;
  b1 <- b1_mu + sig_b1*b1_raw;

  for(n in 1:N){
    mu[n] <- inv_logit(gint[n] + a[yid[n]]  + b1[yid[n]]*X[n] + crowdEff[n]);
  }
  
}
model{
  // Priors
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,4);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);

  // Likelihood
  Y ~ binomial(1,mu);

}
generated quantities {
  vector[Nhold] crowdhat;
  vector[nyrshold] a_out;
  vector[nyrshold] b1_out;
  real muhat[Nhold];
  int<lower=0,upper=1> y_hat[Nhold];          // pointwise predictions  
  vector[Nhold] log_lik;                      // vector for computing log pointwise predictive density  
  vector[Nhold] Whold_intra;          
  
  // 1. Holdout data predictions 
  Whold_intra <- Whold[, spp];
  
  crowdhat <- Whold_intra*w;
  
  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:Nhold){
    muhat[n] <- inv_logit(gint[gidhold[n]] + a_out[yidhold[n] - nyrs] + b1_out[yidhold[n] - nyrs]*Xhold[n] + crowdhat[n]);
    y_hat[n] <- bernoulli_rng(muhat[n]);
    log_lik[n] <- bernoulli_log(Yhold[n], muhat[n]);
  }
  
}
