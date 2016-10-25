// Single-species model for recruitment: includes intraspecific effects 
data{
  int<lower=0> N;                 // observations
  int<lower=0> Y[N];              // observation vector
  int<lower=0> nyrs;              // years
  int<lower=0> yid[N];            // year id
  int<lower=0> G;                 // groups
  int<lower=0> Nspp;              // number of species 
  int<lower=0> spp;               // focal species id
  matrix[N, Nspp] parents1;       // parents in plot
  matrix[N, Nspp] parents2;       // parents in group
  matrix[N, G] gm;
  int<lower=0> Covs;              // climate covariates
  matrix[N,Covs] C;               // climate matrix
  real tau_beta;                  // prior standard deviation
  
  // for out of sample prediction
  int<lower=0> Nhold;
  int<lower=0> nyrshold;              // years out
  int<lower=0> yidhold[Nhold];        // year out id
  int<lower=0> Yhold[Nhold];          // observation vector
  matrix[Nhold, Nspp] parents1hold;   // hold out parents in plot
  matrix[Nhold, Nspp] parents2hold;   // hold out parents in group
  matrix[Nhold, G] gmhold;
  matrix[Nhold,Covs] Chold;           // climate matrix


  // for year effect estimation using entire dataset
  int<lower=0> N2;                    // all observations
  int<lower=0> Y2[N2];                // observation vector
  int<lower=0> nyrs2;                  // all years
  int<lower=0> yid2[N2];              // year id
  matrix[N2, Nspp] parents12;        // parents in plot
  matrix[N2, Nspp] parents22;        // parents in group
  matrix[N2, G] gm2;

}parameters{
  vector[nyrs] a_raw;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  vector[Nspp] w;
  vector[G] bg; 
  vector[Covs] b2; 
  
  // for year effects model
  vector[nyrs2] a_raw2;
  vector[Nspp] w2;
  real<lower=1e-7> sig_a2;
  real<lower=0.0001> theta2;
  real<lower=0, upper=1> u2;
  vector[G] bg2;
}
transformed parameters{
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] q;
  vector[N] coverEff;
  vector[N] gint; 
  vector[nyrs] a; 
  vector[N] climEff; 
  
  // for year effects model
  vector[N2] mu2;
  matrix[N2, Nspp] trueP12;
  matrix[N2, Nspp] trueP22;
  vector[N2] lambda2;
  vector[N2] q2;
  vector[N2] coverEff2;
  vector[N2] gint2;
  vector[nyrs2] a2;

  // for training data model 
  climEff <- C*b2;
  trueP1 <- parents1*u + parents2*(1-u);
  
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  gint     <- gm*bg;
  coverEff <- trueP2*w;
  a  <- 0 + a_raw*sig_a; 

  for(n in 1:N){
    mu[n] <- exp(gint[n] + a[yid[n]] + coverEff[n] + climEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  // elementwise multiplication  
    q[n] <- fmax(lambda[n]*theta, 1e-9);
  }
  
  // for year effects model

  trueP12 <- parents12*u2 + parents22*(1-u2);

  for(n in 1:N2)
    for( j in 1:Nspp)
      trueP22[n, j] <- sqrt(trueP12[n, j]);

  gint2     <- gm2*bg2;
  coverEff2 <- trueP22*w2;
  a2        <- 0 + a_raw2*sig_a2;

  for(n in 1:N2){
    mu2[n] <- exp(gint2[n] + a2[yid2[n]] + coverEff2[n]);
    lambda2[n] <- trueP12[n, spp]*mu2[n];
    q2[n] <- fmax(lambda2[n]*theta2, 1e-9);
  }

}
model{
   // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,5);
  sig_a ~ cauchy(0,5);
  w ~ normal(0, 5);
  a_raw ~ normal(0, 1);
  bg ~ normal(0, 10);
  b2 ~ normal(0, tau_beta);
  
  // Likelihood
  Y ~ neg_binomial_2(q, theta);
  
  // For year effects model
  u2 ~ uniform(0,1);
  theta2 ~ uniform(0,5);
  sig_a2 ~ cauchy(0,5);
  a_raw2 ~ normal(0, 1);
  w2 ~ normal(0, 5);

  // Likelihood
  Y2 ~ neg_binomial_2(q2, theta2);
}
generated quantities{
  vector[nyrshold] a_pred;
  vector[Nhold] coverEff_pred;
  matrix[Nhold, Nspp] trueP1_pred;
  matrix[Nhold, Nspp] trueP2_pred;
  vector[Nhold] mu_pred;
  vector[Nhold] lambda_pred;
  vector[Nhold] q_pred;
  vector[Nhold] log_lik; // vector for computing log pointwise predictive density
  //int<lower=0> y_hat[Nhold]; // pointwise predictions  
  vector[Nhold] gint_pred; 
  vector[Nhold] climEff_pred; 
  
  //for year predictions from year effects model
  vector[Nhold] q_pred2;
  vector[Nhold] mu_pred2;
  vector[Nhold] lambda_pred2;
  //int<lower=0> y_hat2[Nhold]; // pointwise predictions
  vector[Nhold] log_lik2; // vector for computing log pointwise predictive density

  // 1. Holdout data predictions 

  climEff_pred <- Chold*b2;
  gint_pred   <- gmhold*bg;
  trueP1_pred <- parents1hold*u + parents2hold*(1-u);

  for(n in 1:Nhold)
    for(j in 1:Nspp)
      trueP2_pred[n, j] <- sqrt(trueP1_pred[n, j]);

  coverEff_pred <- trueP2_pred*w;

  for( i in 1:nyrshold)
    a_pred[i] <- normal_rng(0, sig_a); // draw random year intercept

  for(n in 1:Nhold){
    mu_pred[n] <- exp(gint_pred[n] + a_pred[yidhold[n] - nyrs ] + coverEff_pred[n] + climEff_pred[n]);
    lambda_pred[n] <- trueP1_pred[n, spp]*mu_pred[n];
    q_pred[n] <- fmax( lambda_pred[n]*theta, 1e-9);
  }

  for(n in 1:Nhold){
    //y_hat[n] <- neg_binomial_2_rng(q_pred[n],  theta);
    log_lik[n] <- neg_binomial_2_log(Yhold[n], q_pred[n], theta);
  }

  // 2. Predictions for holdout data with KNOWN year effects.  
  //    Simulate predictions as if year effects in the out of sample data are known. 
  
  for(n in 1:Nhold){
    mu_pred2[n] <- exp(gint_pred[n] + a2[yidhold[n]] + coverEff_pred[n]);
    lambda_pred2[n] <- trueP1_pred[n, spp]*mu_pred2[n];  // elementwise multiplication 
    q_pred2[n] <- fmax( lambda_pred2[n]*theta, 1e-9);
  }
  
  for(n in 1:Nhold){
    //y_hat2[n] <- neg_binomial_2_rng(q_pred2[n],  theta);
    log_lik2[n] <- neg_binomial_2_log(Yhold[n], q_pred2[n], theta);
  } 
  
}
