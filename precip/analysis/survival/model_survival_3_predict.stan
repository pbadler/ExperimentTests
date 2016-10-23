data{
  // training datalist, historical observations 
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
  int<lower=0> Covs;          // climate effects 
  matrix[N, Covs] C;          // climate matrix 
  real <lower=0> tau_beta;    // prior sd for climate effects 
  
  // full datalist, all observations for year effects  
  int<lower=0> N2;              // observations
  int<lower=0,upper=1> Y2[N2];  // observation vector
  vector[N2] X2;                // size vector
  matrix[N2, G] gm2;            // group dummy variable matrix 
  int<lower=0> nyrs2;           // years
  int<lower=0> yid2[N2];        // year id
  matrix[N2,Wcovs] W2;          // crowding matrix

  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  int<lower=0,upper=1> Yhold[Nhold];
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  matrix[Nhold, G] gmhold;          // group dummy variable matrix 
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
  matrix[Nhold, Covs] Chold;        // climate matrix 

}
parameters{
  // for training data model  
  vector[G] bg;                     // varying group effects with first group as intercept 
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[Wcovs] w;
  vector[Covs]  b2; 
  
  // for year effects model 
  vector[G] bg2;                     // varying group effects with first group as intercept 
  vector[nyrs2] a_raw2;
  real b1_mu2;
  vector[nyrs2] b1_raw2;
  real<lower=0> sig_a2;
  real<lower=0> sig_b12;
  vector[Wcovs] w2;
}
transformed parameters{
  // for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff; 
  
  // for year effects model 
  vector[nyrs2] a2;
  vector[nyrs2] b12;
  vector[N2] gint2; 
  real mu2[N2];
  vector[N2] crowdEff2;

  // for training data model -----------------------------------
  gint <- gm*bg;
  crowdEff <- W*w;
  climEff  <- C*b2;
  
  // reparamaterize the year effects parameters  
  a <- 0 + sig_a*a_raw;
  b1 <- b1_mu + sig_b1*b1_raw;

  for(n in 1:N){
    mu[n] <- inv_logit(gint[n] + a[yid[n]]  + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n]);
  }

  // for year effects model -----------------------------------
  gint2 <- gm2*bg2;
  crowdEff2 <- W2*w2;
  
  // reparamaterize the year effects parameters  
  a2 <- 0 + sig_a2*a_raw2;
  b12 <- b1_mu2 + sig_b12*b1_raw2;

  for(n in 1:N2){
    mu2[n] <- inv_logit(gint2[n] + a2[yid2[n]]  + b12[yid2[n]]*X2[n] + crowdEff2[n]);
  }
  
}
model{
  // for training data model 
  // Priors
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,4);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);
  b2 ~ normal(0,tau_beta); 
    
  // Likelihood
  Y ~ binomial(1,mu);
  
  // for year effects model 
  // Priors
  bg2 ~ normal(0,10);
  b1_mu2 ~ normal(0,10);
  sig_a2 ~ cauchy(0,5);
  sig_b12 ~ cauchy(0,4);
  a_raw2 ~ normal(0,1);
  b1_raw2 ~ normal(0,1);
  w2 ~ normal(0,10);

  // Likelihood
  Y2 ~ binomial(1,mu2);

}
generated quantities {
  // hold out predictions 
  vector[nyrshold] a_out;
  vector[nyrshold] b1_out;
  vector[Nhold] gint_out;
  real muhat[Nhold];
  int<lower=0,upper=1> y_hat[Nhold];          // pointwise predictions  
  vector[Nhold] log_lik;                      // vector for computing log pointwise predictive density  
  vector[Nhold] crowdhat;
  vector[Nhold] climhat;
  
  // for predictions with year effects 
  real muhat2[Nhold];
  vector[Nhold] gint_out2;
  int<lower=0,upper=1> y_hat2[Nhold];          // pointwise predictions  
  vector[Nhold] log_lik2;                      // vector for computing log pointwise predictive density  
  
  // 1. Holdout data predictions 
  gint_out  <- gmhold*bg;
  crowdhat <- Whold*w;
  climhat  <- Chold*b2; 
  
  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:Nhold){
    muhat[n] <- inv_logit(gint_out[n] + a_out[yidhold[n] - nyrs] + b1_out[yidhold[n] - nyrs]*Xhold[n] + crowdhat[n] + climhat[n]);
    y_hat[n] <- bernoulli_rng(muhat[n]);
    log_lik[n] <- bernoulli_log(Yhold[n], muhat[n]);
  }
  
  // 2. Holdout data predictions with known year effects 
  gint_out2 <- gmhold*bg;
  
  for(n in 1:Nhold){
    muhat2[n] <- inv_logit(gint_out2[n] + a2[yidhold[n]] + b12[yidhold[n]]*Xhold[n] + crowdhat[n]);
    y_hat2[n] <- bernoulli_rng(muhat2[n]);
    log_lik2[n] <- bernoulli_log(Yhold[n], muhat2[n]);
  }
  
}
