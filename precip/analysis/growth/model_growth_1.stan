data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
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
  real<lower=0> tau_beta;    // prior sd for climate effects
  
  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  vector[Nhold] Yhold;
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
  real sigma; 
}
transformed parameters{
  // for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff; 

  // for training data model -----------------------------------
  gint <- gm*bg;
  crowdEff <- W*w;
  climEff  <- C*b2;
  
  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n];
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
  sigma ~ cauchy(0, 5);
  w ~ normal(0,10);
  b2 ~ normal(0,tau_beta); 
    
  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {
  // hold out predictions
  vector[N] log_lik;                          // for fitted data
  vector[Nhold] log_lik2;                     // for heldout data 

  vector[nyrshold] a_out;
  vector[nyrshold] b1_out;
  real muhat[Nhold];
  vector[Nhold] gint_out;
  vector[Nhold] crowdhat;
  vector[Nhold] climhat;

  # fitted data log_lik 
  for(n in 1:N){
      log_lik[n] <- normal_log(Y[n], mu[n], sigma);
  }
    
  // 1. Holdout data predictions
  gint_out  <- gmhold*bg;
  crowdhat  <- Whold*w;
  climhat   <- Chold*b2;

  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept
    b1_out[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect
  }

  for(n in 1:Nhold){
      muhat[n]    <- gint_out[n] + a_out[yidhold[n]-nyrs] + b1_out[yidhold[n]-nyrs]*Xhold[n] + crowdhat[n] + climhat[n];
      log_lik2[n]  <- normal_log(Yhold[n], muhat[n], sigma);
  }

}
