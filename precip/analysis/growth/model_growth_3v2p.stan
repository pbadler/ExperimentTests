data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  matrix[N, 2] X;                // size vector
  int<lower=0> G;             // groups
  matrix[N, G] gm;            // group dummy variable matrix 
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  int<lower=0> Wcovs;         // number of crowding effects 
  matrix[N,Wcovs] W;          // crowding matrix
  int<lower=0> spp;           // focal species 
  int<lower=0> Covs;          // climate effects 
  matrix[N, Covs] C;          // climate matrix 
  real<lower=0> tau_beta;    // prior sd for climate effects
  
  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  vector[Nhold] Yhold;
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  matrix[Nhold, G] gmhold;          // group dummy variable matrix
  matrix[Nhold, 2] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
  matrix[Nhold, Covs] Chold;        // climate matrix
}
parameters{
  // for training data model  
  vector[G] bg;                     // varying group effects with first group as intercept 
  vector[nyrs] a_raw;
  vector[2] b1_mu;
  vector[nyrs] b1_raw1;
  vector[nyrs] b1_raw2;
  real<lower=0> sig_a;
  vector<lower=0>[2] sig_b1;
  vector[Wcovs] w;
  vector[Covs]  b2; 
  real<lower=0> tau; 
  real<lower=0> tauSize; 
}
transformed parameters{
  // for training data model  
  vector[nyrs] a;
  vector[nyrs] b11;
  vector[nyrs] b12;
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff; 
  vector[N] sigma;
  
  // for training data model -----------------------------------
  gint <- gm*bg;
  crowdEff <- W*w;
  climEff  <- C*b2;
  
  b11 <- b1_mu[1] + sig_b1[1]*b1_raw1;
  b12 <- b1_mu[2] + sig_b1[2]*b1_raw2;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b11[yid[n]]*X[n, 1] + b12[yid[n]]*X[n, 2] + crowdEff[n] + climEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
  }
}
model{
  // Priors
  bg ~ normal(0,5);
  b1_mu ~ normal(0,2);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw1 ~ normal(0,1);
  
  tau ~ normal( 0, 1); 
  tauSize ~ normal(0,1);
  w ~ normal(0,tau_beta);
  b2 ~ normal(0,tau_beta); 
    
  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {
  // hold out predictions
  vector[N] log_lik;                          // for fitted data
  vector[Nhold] log_lik2;                     // for heldout data 

  vector[nyrshold] a_out;
  vector[nyrshold] b1_out1;
  vector[nyrshold] b1_out2;

  real muhat[Nhold];
  vector[Nhold] gint_out;
  vector[Nhold] y_hat;                        // pointwise predictions
  vector[Nhold] crowdhat;
  vector[Nhold] climhat;
  vector[Nhold] sigmahat;

  # fitted data log_lik 
  for(n in 1:N){
      log_lik[n] <- normal_log(Y[n], mu[n], sigma[n]);
  }
    
  // 1. Holdout data predictions
  gint_out  <- gmhold*bg;
  crowdhat  <- Whold*w;
  climhat   <- Chold*b2;

  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept
    b1_out1[i] <- normal_rng(b1_mu[1], sig_b1[1]);   //draw random year x size effect
    b1_out2[i] <- normal_rng(b1_mu[2], sig_b1[2]);   //draw random year x size effect
  }

  for(n in 1:Nhold){
      muhat[n]    <- gint_out[n] + a_out[yidhold[n]-nyrs] + b1_out1[yidhold[n]-nyrs]*Xhold[n, 1] + b1_out2[yidhold[n]-nyrs]*Xhold[n, 2] + crowdhat[n] + climhat[n];
      sigmahat[n] <- sqrt((fmax(tau*exp(tauSize*muhat[n]), 0.0000001)));
      y_hat[n]    <- normal_rng(muhat[n], sigmahat[n]);
      log_lik2[n]  <- normal_log(Yhold[n], muhat[n], sigmahat[n]);
  }

}
