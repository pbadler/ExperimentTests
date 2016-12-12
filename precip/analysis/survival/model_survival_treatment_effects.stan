data{
  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  int<lower=0, upper=1> Yhold[Nhold];
  int<lower=0> nyrs;                // training data years 
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  int<lower=0> G;                   // Groups 
  matrix[Nhold, G] gmhold;          // group dummy variable matrix
  int<lower=0> nThold;              // treatment groups
  matrix[Nhold, nThold] tmhold;     // group dummy variable matrix
  int<lower=0> Wcovs;               // number of crowding effects 
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
  
  // holdout datalist for cover predictions 
  int<lower=0> N2;
  int<lower=0> nyrs2;            // years out
  int<lower=0> yid2[N2];      //year out id
  matrix[N2, G] gm2;          // group dummy variable matrix 
  vector[N2] X2;
  matrix[N2,Wcovs] W2;        // crowding matrix for holdout data
  matrix[N2, nThold] tm2;     // treatment matrix 
}
parameters{
  // for training data model  
  real b1_mu;
  vector[nyrshold] a_raw;
  vector[nyrshold] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[Wcovs] w;
  vector[nThold] bt;
  vector[G] bg;
}
transformed parameters{
  // for training data model  
  vector[Nhold] treatEff;
  real mu[Nhold];
  vector[Nhold] crowdEff;
  vector[Nhold] gint;
  vector[nyrshold] a;
  vector[nyrshold] b1;

  // for training data model -----------------------------------
  crowdEff <- Whold*w;
  treatEff <- tmhold*bt;
  gint     <- gmhold*bg;
  
  a <- 0 + sig_a*a_raw; 
  b1 <- b1_mu + sig_b1*b1_raw;
  
  for(n in 1:Nhold){
    mu[n] <- inv_logit(gint[n] + treatEff[n] + a[yidhold[n] - nyrs] + b1[yidhold[n] - nyrs]*Xhold[n] + crowdEff[n]);
  }
  
}
model{
  // Priors
  bt ~ normal(0,10);
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0, 1);
  b1_raw ~ normal(0, 1);
  w ~ normal(0,10);

  // Likelihood
  Yhold ~ bernoulli_log(mu);
}
generated quantities {
  vector[Nhold] log_lik2;                     // for heldout data 
  vector[Nhold] year_effect; 
  
  // for cover predictions 
  vector[nyrs2] a_out2;
  vector[nyrs2] b1_out2;
  vector[N2] gint_out2;
  real muhat2[N2];
  vector[N2] crowdhat2;
  vector[N2] treatEff2;

  for(n in 1:Nhold)
    year_effect[n] <-  a[yidhold[n] - nyrs] + Xhold[n]*(b1[yidhold[n] - nyrs] - b1_mu);
  
  for(n in 1:Nhold){
    log_lik2[n] <- bernoulli_log(Yhold[n], mu[n]); 
  }
  
  // 2. all data for cover predictions 
  gint_out2  <- gm2*bg;
  crowdhat2 <- W2*w;
  treatEff2 <- tm2*bt; 
  
  for( i in 1:nyrs2){
    a_out2[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out2[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:N2){
    muhat2[n] <- inv_logit(gint_out2[n] + a_out2[yid2[n]] + b1_out2[yid2[n]]*X2[n] + crowdhat2[n] + treatEff2[n]);
  }
  
  
}

