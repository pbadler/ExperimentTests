data{
  int<lower=0> N;
  int<lower=0, upper=1> Y[N];
  int<lower=0> nyrs;                // years out
  int<lower=0> yid[N];              //year out id
  int<lower=0> G;                   // Groups 
  matrix[N, G] gm;                  // group dummy variable matrix
  int<lower=0> Wcovs;               // number of crowding effects 
  vector[N] X;
  matrix[N,Wcovs] W;                // crowding matrix for holdout data
  
  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  int<lower=0,upper=1> Yhold[Nhold];
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  matrix[Nhold, G] gmhold;          // group dummy variable matrix 
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
  
  // holdout datalist for cover predictions 
  int<lower=0> N2;
  int<lower=0> nyrs2;            // years out
  int<lower=0> yid2[N2];      //year out id
  matrix[N2, G] gm2;          // group dummy variable matrix 
  vector[N2] X2;
  matrix[N2,Wcovs] W2;        // crowding matrix for holdout data
}
parameters{
  // for training data model  
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[Wcovs] w;
  vector[G] bg;
}
transformed parameters{
  // for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  real mu[N];
  vector[N] crowdEff;
  vector[N] gint;

  // for training data model -----------------------------------
  crowdEff <- W*w;
  gint     <- gm*bg;
  
  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- inv_logit(gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n]);
  }
  
}
model{
  // Priors
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);

  // Likelihood
  Y ~ bernoulli_log(mu);
}
generated quantities {
  // hold out predictions
  vector[N] log_lik;                          // for fitted data
  vector[Nhold] log_lik2;                     // for heldout data 
  
  // hold out predictions 
  vector[nyrshold] a_out;
  vector[nyrshold] b1_out;
  vector[Nhold] gint_out;
  real muhat[Nhold];
  vector[Nhold] crowdhat;

  // for cover predictions 
  vector[nyrs2] a_out2;
  vector[nyrs2] b1_out2;
  vector[N2] gint_out2;
  real muhat2[N2];
  vector[N2] crowdhat2;

  for(n in 1:N){
    log_lik[n] <- bernoulli_log(Y[n], mu[n]); 
  }
  
  // 1. Holdout data predictions 
  gint_out  <- gmhold*bg;
  crowdhat <- Whold*w;

  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:Nhold){
    muhat[n] <- inv_logit(gint_out[n] + a_out[yidhold[n]-nyrs] + b1_out[yidhold[n]-nyrs]*Xhold[n] + crowdhat[n]);
    log_lik2[n] <- bernoulli_log(Yhold[n], muhat[n]);
  }
  
  // 2. all data for cover predictions 
  gint_out2  <- gm2*bg;
  crowdhat2 <- W2*w;

  for( i in 1:nyrs2){
    a_out2[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out2[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:N2){
    muhat2[n] <- inv_logit(gint_out2[n] + a_out2[yid2[n]] + b1_out2[yid2[n]]*X2[n] + crowdhat2[n]);
  }
  
  
}
