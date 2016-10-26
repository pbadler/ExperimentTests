// Single-species model for recruitment: includes intraspecific effects 
data{
  int<lower=0> N;             // observations
  int<lower=0> Y[N];          // observation vector
  int<lower=0> nyrs;           // years
  int<lower=0> yid[N];        // year id
  int<lower=0> G;             // groups
  int<lower=0> Nspp;          // number of species 
  int<lower=0> spp;           // focal species id
  matrix[N, Nspp] parents1;   // parents in plot
  matrix[N, Nspp] parents2;   // parents in group
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
  real w;
  vector[G] bg; 
  vector[Covs] b2; 
  
  // for year effects model  
  vector [nyrs2] a_raw2;
  real w2;
  real<lower=1e-7> sig_a2;
  real<lower=0> theta2;
  real<lower=0, upper=1> u2;
  vector[G] bg2; 

}
transformed parameters{
  vector[N] mu;
  vector[N] trueP1;
  vector[N] trueP2;
  vector[N] lambda;
  vector[N] coverEff;
  vector[N] p1; 
  vector[N] p2;
  vector[N] gint; 
  vector[nyrs] a; 
  vector[N] climEff; 
  
  // for year effects model 
  vector[N2] mu2; 
  vector[N2] trueP1_2;
  vector[N2] trueP2_2;
  vector[N2] lambda2;
  vector[N2] coverEff2;
  vector[N2] p1_2; 
  vector[N2] p2_2;
  vector[N2] gint2; 
  vector[nyrs2] a2; 
  
  for( n in 1:N){ 
    p1[n] <- parents1[n, spp];
    p2[n] <- parents2[n, spp];
  }
  
  trueP1 <- p1*u + p2*(1-u);

  for(n in 1:N)
      trueP2[n] <- sqrt(trueP1[n]);
      
  climEff   <- C*b2;
  gint      <- gm*bg;
  coverEff  <- trueP2*w;
  a         <- 0 + a_raw*sig_a; 
  
  for(n in 1:N){
    mu[n] <- exp(gint[n]  + a[yid[n]]  + coverEff[n] + climEff[n]);
    lambda[n] <- trueP1[n]*mu[n];  
  } 
  
  // for year effects model 
  for( n in 1:N2){ 
    p1_2[n] <- parents12[n, spp];
    p2_2[n] <- parents22[n, spp];
  }
  
  trueP1_2 <- p1_2*u2 + p2_2*(1-u2);

  for(n in 1:N2)
      trueP2_2[n] <- sqrt(trueP1_2[n]);
  
  gint2     <- gm2*bg2;
  coverEff2 <- trueP2_2*w2;
  a2        <- 0 + a_raw2*sig_a2; 

  for(n in 1:N2){
    mu2[n] <- exp(gint2[n]  + a2[yid2[n]]  + coverEff2[n]);
    lambda2[n] <- trueP1_2[n]*mu2[n];  // elementwise multiplication  
  }

}
model{
   // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,2);
  sig_a ~ cauchy(0,2);
  a_raw ~ normal(0, 1);
  bg ~ normal(0, 10);
  w ~ normal(0, 5);
  b2 ~ normal(0, tau_beta); 
  
  // Likelihood
  Y ~ neg_binomial_2(lambda, theta);
  
  // For year effects model 
  u2 ~ uniform(0,1);
  theta2 ~ cauchy(0,2);
  sig_a2 ~ cauchy(0,2);
  a_raw2 ~ normal(0, 1);
  bg2 ~ normal(0, 10);
  w2 ~ normal(0, 5);

  // Likelihood
  Y2 ~ neg_binomial_2(lambda2, theta2);
}
generated quantities{
  vector[nyrshold] a_pred;
  vector[Nhold] coverEff_pred;
  vector[Nhold] trueP1_pred;
  vector[Nhold] trueP2_pred;
  vector[Nhold] mu_pred;
  vector[Nhold] lambda_pred;
  vector[Nhold] log_lik; // vector for computing log pointwise predictive density
  int<lower=0> y_hat[Nhold]; // pointwise predictions  
  vector[Nhold] p1_pred; 
  vector[Nhold] p2_pred;
  vector[Nhold] gint_pred; 
  vector[Nhold] climEff_pred; 
  
  // for year predictions from year effects model 
  vector[Nhold] mu_pred2;
  vector[Nhold] lambda_pred2;
  int<lower=0> y_hat2[Nhold]; // pointwise predictions  
  vector[Nhold] log_lik2; // vector for computing log pointwise predictive density  

  // 1. Holdout data predictions 
  for( n in 1:Nhold){ 
    p1_pred[n] <- parents1hold[n, spp];
    p2_pred[n] <- parents2hold[n, spp];
  }
  
  climEff_pred<- Chold*b2;
  gint_pred   <- gmhold*bg;
  trueP1_pred <- p1_pred*u + p2_pred*(1-u);

  for(n in 1:Nhold)
      trueP2_pred[n] <- sqrt(trueP1_pred[n]);
  
  coverEff_pred <- trueP2_pred*w;

  for( i in 1:nyrshold)
    a_pred[i] <- normal_rng(0, sig_a); // draw random year intercept 

  for(n in 1:Nhold){
    mu_pred[n] <- exp(gint_pred[n] + a_pred[yidhold[n] - nyrs ] + coverEff_pred[n] + climEff_pred[n]);
    lambda_pred[n] <- trueP1_pred[n]*mu_pred[n];
  }
  
  for(n in 1:Nhold){
    y_hat[n] <- neg_binomial_2_rng(lambda_pred[n],  theta);
    log_lik[n] <- neg_binomial_2_log(Yhold[n], lambda_pred[n], theta);
  }
  
  // 2. Predictions for holdout data with KNOWN year effects.  
  //    Simulate predictions as if year effects in the out of sample data are known. 
  
  for(n in 1:Nhold){
    mu_pred2[n] <- exp(gint_pred[n] + a2[yidhold[n]] + coverEff_pred[n]);
    lambda_pred2[n] <- trueP1_pred[n]*mu_pred2[n];  // elementwise multiplication 
  }
  
  for(n in 1:Nhold){
    y_hat2[n] <- neg_binomial_2_rng(lambda_pred2[n],  theta);
    log_lik2[n] <- neg_binomial_2_log(Yhold[n], lambda_pred2[n], theta);
  }  
}
