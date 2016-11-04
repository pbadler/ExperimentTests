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
  
  // // holdout datalist, modern observations 
  int<lower=0> Nhold;
  vector[Nhold] Yhold;
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  matrix[Nhold, G] gmhold;          // group dummy variable matrix
  vector[Nhold] Xhold;
  matrix[Nhold,Wcovs] Whold;        // crowding matrix for holdout data
  matrix[Nhold, Covs] Chold;        // climate matrix

  // survival data for cover predictions
  int<lower=0> N3;              // observations
  vector[N3] X3;                // size vector
  matrix[N3, G] gm3;            // group dummy variable matrix
  int<lower=0> nyrs3;           // years
  int<lower=0> yid3[N3];        // year id
  matrix[N3, Wcovs] W3;         // crowding matrix
  matrix[N3, Covs] C3;          // climate matrix
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
  real<lower=0> tau; 
  real<lower=0> tauSize; 
}
transformed parameters{
  // for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff; 
  vector[N] sigma;
  
  // for training data model -----------------------------------
  gint <- gm*bg;
  crowdEff <- W*w;
  climEff  <- C*b2;
  
  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
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
  tau ~ normal( 0, 1); 
  tauSize ~ normal(0,1);
  w ~ normal(0,tau_beta);
  b2 ~ normal(0,tau_beta); 
    
  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {
  // hold out predictions
  vector[nyrshold] a_out;
  vector[nyrshold] b1_out;
  real muhat[Nhold];
  vector[Nhold] gint_out;
  vector[Nhold] y_hat;                        // pointwise predictions
  vector[Nhold] log_lik;                      // vector for computing log pointwise predictive density
  vector[Nhold] crowdhat;
  vector[Nhold] climhat;
  vector[Nhold] sigmahat;

  // for cover predictions
  vector[N3] crowdhat3;
  vector[N3] climhat3;
  vector[N3] gint3;
  real muhat3[N3];
  vector[N3] y_hat3;                          // pointwise predictions
  vector[N3] sigmahat3;

  // 1. Holdout data predictions
  gint_out  <- gmhold*bg;
  crowdhat  <- Whold*w;
  climhat   <- Chold*b2;

  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept
    b1_out[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect
  }

  for(n in 1:Nhold){
      muhat[n]    <- gint_out[n] + a[yidhold[n]-nyrs] + b1[yidhold[n]-nyrs]*Xhold[n] + crowdhat[n] + climhat[n];
      sigmahat[n] <- sqrt((fmax(tau*exp(tauSize*muhat[n]), 0.0000001)));
      y_hat[n]    <- normal_rng(muhat[n], sigmahat[n]);
      log_lik[n]  <- normal_log(Yhold[n], muhat[n], sigmahat[n]);
  }

  // 3. Holdout survival data for predicting cover
  gint3     <- gm3*bg;
  crowdhat3 <- W3*w;
  climhat3  <- C3*b2;

  for(n in 1:N3){
      muhat3[n] <- gint3[n] + a_out[yid3[n]-nyrs] + b1_out[yid3[n]-nyrs]*X3[n] + crowdhat3[n] + climhat3[n];
      sigmahat3[n] <- sqrt((fmax(tau*exp(tauSize*muhat3[n]), 0.0000001)));
      y_hat3[n] <- normal_rng(muhat3[n], sigmahat3[n]);
  }
}
