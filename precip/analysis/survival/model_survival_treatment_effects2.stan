data{
  // holdout datalist, modern observations 
  int<lower=0> N2;
  int<lower=0, upper=1> Y2[N2];
  int<lower=0> nyrs2;            // years out
  int<lower=0> yid2[N2];      //year out id
  int<lower=0> G;                   // Groups 
  matrix[N2, G] gm2;          // group dummy variable matrix
  int<lower=0> nT;              // treatment groups
  matrix[N2, nT] tm2;     // group dummy variable matrix
  int<lower=0> Wcovs;               // number of crowding effects 
  vector[N2] X2;
  matrix[N2,Wcovs] W2;        // crowding matrix for holdout data
  
}
parameters{
  // for training data model  
  real b1_mu;
  vector[nyrs2] a_raw;
  vector[nyrs2] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[Wcovs] w;
  vector[nT] bt;
  vector[G] bg;
}
transformed parameters{
  // for training data model  
  vector[N2] treatEff;
  real mu[N2];
  vector[N2] crowdEff;
  vector[N2] gint;
  vector[nyrs2] a;
  vector[nyrs2] b1;

  // for training data model -----------------------------------
  crowdEff <- W2*w;
  treatEff <- tm2*bt;
  gint     <- gm2*bg;
  
  a <- 0 + sig_a*a_raw; 
  b1 <- b1_mu + sig_b1*b1_raw;
  
  for(n in 1:N2){
    mu[n] <- inv_logit(gint[n] + treatEff[n] + a[yid2[n]] + b1[yid2[n]]*X2[n] + crowdEff[n]);
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
  Y2 ~ bernoulli_log(mu);
}
generated quantities {
  vector[N2] log_lik2;                     // for heldout data 
  vector[N2] year_effect; 
  
  // for cover predictions 
  vector[nyrs2] a_out2;
  vector[nyrs2] b1_out2;
  real muhat2[N2];

  for(n in 1:N2)
    year_effect[n] <-  a[yid2[n]] + X2[n]*(b1[yid2[n]] - b1_mu);
  
  for(n in 1:N2){
    log_lik2[n] <- bernoulli_log(Y2[n], mu[n]); 
  }
  
  // 2. all data for cover predictions 
  for( i in 1:nyrs2){
    a_out2[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out2[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:N2){
    muhat2[n] <- inv_logit(gint[n] + a_out2[yid2[n]] + b1_out2[yid2[n]]*X2[n] + crowdEff[n] + treatEff[n]);
  }
  
  
}

