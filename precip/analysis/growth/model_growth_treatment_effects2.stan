data{
  // holdout datalist, modern observations 
  int<lower=0> N2;
  vector[N2] Y2;
  int<lower=0> nyrs2;                // training data years 
  int<lower=0> yid2[N2];      // year out id
  int<lower=0> G;                   // Groups 
  matrix[N2, G] gm2;          // group dummy variable matrix
  int<lower=0> nT;              // number of treatments 
  matrix[N2, nT] tm2;     // treatment dummy matrix  
  int<lower=0> Wcovs;               // number of crowding effects 
  vector[N2] X2;
  matrix[N2,Wcovs] W2;        // crowding matrix for holdout data
  
  // holdout datalist for cover 
  int<lower=0> N3;
  vector[N3] Y3;
  int<lower=0> nyrs3;         // years out
  int<lower=0> yid3[N3];      // year out id
  matrix[N3, G] gm3;          // group dummy variable matrix
  vector[N3] X3;
  matrix[N3,Wcovs] W3;        // crowding matrix for holdout data
  matrix[N3, nT] tm3;     // treatment dummy matrix  

}
parameters{
  // for training data model  
  real b1_mu;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[nyrs2] a_raw;
  vector[nyrs2] b1_raw;
  real<lower=0> sigma;
  vector[Wcovs] w;
  vector[nT] bt;
  vector[G] bg;
}
transformed parameters{
  // for training data model  
  vector[N2] treatEff;
  vector[N2] gint;
  real mu[N2];
  vector[N2] crowdEff;
  vector[nyrs2] a; 
  vector[nyrs2] b1; 
  
  // for training data model -----------------------------------
  crowdEff <- W2*w;
  treatEff <- tm2*bt;
  gint <- gm2*bg; 
  
  a <- 0 + sig_a*a_raw; 
  b1 <- b1_mu + sig_b1*b1_raw;
  
  for(n in 1:N2){
    mu[n] <- gint[n] + a[yid2[n]] + b1[yid2[n]]*X2[n]+ treatEff[n]  + crowdEff[n];
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
  Y2 ~ normal(mu, sigma);
}
generated quantities {
  vector[N2] log_lik2;                     // for heldout data 
  vector[N2] year_effect; 
  
  vector[nyrs3] a_out3;
  vector[nyrs3] b1_out3;
  real muhat3[N3];
  vector[N3] gint_out3;
  vector[N3] crowdhat3;
  vector[N3] treatEff3;
  
  for(n in 1:N2)
    year_effect[n] <-  a[yid2[n]] + X2[n]*(b1[yid2[n]] - b1_mu);
  
  # fitted data log_lik 
  for(n in 1:N2){
      log_lik2[n] <- normal_log(Y2[n], mu[n], sigma);
  }
  
  // 2. cover data predictions 
  gint_out3  <- gm3*bg;
  crowdhat3  <- W3*w;
  treatEff3   <- tm3*bt;

  for( i in 1:nyrs3){
    a_out3[i] <- normal_rng(0, sig_a);         // draw random year intercept
    b1_out3[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect
  }

  for(n in 1:N3){
      muhat3[n]    <- gint_out3[n] + a_out3[yid3[n]] + b1_out3[yid3[n]]*X3[n] + crowdhat3[n] + treatEff3[n];
  }
  
}
