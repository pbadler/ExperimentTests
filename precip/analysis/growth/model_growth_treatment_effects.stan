data{
  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  vector[Nhold] Yhold;
  int<lower=0> nyrs;                // training data years 
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      // year out id
  int<lower=0> G;                   // Groups 
  matrix[Nhold, G] gmhold;          // group dummy variable matrix
  int<lower=0> nT;              // number of treatments 
  matrix[Nhold, nT] tmhold;     // treatment dummy matrix  
  int<lower=0> Wcovs;               // number of crowding effects 
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
  
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
  vector[nyrshold] a_raw;
  vector[nyrshold] b1_raw;
  real<lower=0> sigma;
  vector[Wcovs] w;
  vector[nT] bt;
  vector[G] bg;
}
transformed parameters{
  // for training data model  
  vector[Nhold] treatEff;
  vector[Nhold] gint;
  real mu[Nhold];
  vector[Nhold] crowdEff;
  vector[nyrshold] a; 
  vector[nyrshold] b1; 
  
  // for training data model -----------------------------------
  crowdEff <- Whold*w;
  treatEff <- tmhold*bt;
  gint <- gmhold*bg; 
  
  a <- 0 + sig_a*a_raw; 
  b1 <- b1_mu + sig_b1*b1_raw;
  
  for(n in 1:Nhold){
    mu[n] <- gint[n] + a[yidhold[n] - nyrs ] + b1[yidhold[n] - nyrs]*Xhold[n]+ treatEff[n]  + crowdEff[n];
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
  vector[Nhold] year_effect; 
  
  vector[nyrs3] a_out3;
  vector[nyrs3] b1_out3;
  real muhat3[N3];
  vector[N3] gint_out3;
  vector[N3] crowdhat3;
  vector[N3] treatEff3;
  
  for(n in 1:Nhold)
    year_effect[n] <-  a[yidhold[n] - nyrs] + Xhold[n]*(b1[yidhold[n] - nyrs] - b1_mu);
  
  # fitted data log_lik 
  for(n in 1:Nhold){
      log_lik2[n] <- normal_log(Yhold[n], mu[n], sigma);
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
