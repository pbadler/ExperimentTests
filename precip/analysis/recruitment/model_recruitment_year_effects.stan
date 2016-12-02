data{
  int<lower=0> N;                 // observations
  int<lower=0> Y[N];              // observation vector
  int<lower=0> nyrs;              // years
  int<lower=0> yid[N];            // year id
  int<lower=0> G;                 // Groups 
  matrix[N, G] gm;                // group dummy variable matrix
  int<lower=0> Nspp;              // number of species 
  int<lower=0> spp;               // focal species id
  matrix[N, Nspp] parents1;       // parents in plot
  matrix[N, Nspp] parents2;       // parents in group
  
  // for out of sample prediction
  int<lower=0> Nhold;
  int<lower=0> nyrshold;              // years out
  int<lower=0> yidhold[Nhold];        // year out id
  int<lower=0> Yhold[Nhold];          // observation vector
  matrix[Nhold, Nspp] parents1hold;   // hold out parents in plot
  matrix[Nhold, Nspp] parents2hold;   // hold out parents in group
  matrix[Nhold, G] gmhold;
  
}parameters{
  vector[nyrs] a_raw;
  vector[Nspp] w;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  vector[G] bg;
}
transformed parameters{
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] coverEff;
  vector[N] treatEff;
  vector[nyrs] a; 
  vector[N] gint; 
  
  trueP1 <- parents1*u + parents2*(1-u);
  
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  gint     <- gm*bg;
  coverEff <- trueP2*w;
  a  <- 0 + a_raw*sig_a; 

  for(n in 1:N){
    mu[n] <- exp(gint[n] + a[yid[n]] + coverEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  
  }
}

model{
  // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,2);
  sig_a ~ cauchy(0,2);
  a_raw ~ normal(0, 1);
  bg ~ normal(0, 10);
  w ~ normal(0, 10);
  
  // Likelihood
  Y ~ neg_binomial_2(lambda, theta);
}
generated quantities{
  vector[N] log_lik; 
  vector[Nhold] log_lik2; 
  vector[nyrshold] a_pred;
  vector[Nhold] mu_pred;
  vector[Nhold] lambda_pred;

  // for prediction 
  vector[Nhold] coverEff_pred;
  matrix[Nhold, Nspp] trueP1_pred;
  matrix[Nhold, Nspp] trueP2_pred;
  vector[Nhold] gint_out; 

  for(n in 1:N){ 
    log_lik[n] <- neg_binomial_2_log(Y[n], lambda[n], theta); 
  }
  
  // 1. Holdout data predictions 
  gint_out   <- gmhold*bg;
  trueP1_pred <- parents1hold*u + parents2hold*(1-u);

  for(n in 1:Nhold)
    for(j in 1:Nspp)
      trueP2_pred[n, j] <- sqrt(trueP1_pred[n, j]);
  
  coverEff_pred <- trueP2_pred*w;

  for( i in 1:nyrshold)
    a_pred[i] <- normal_rng(0, sig_a); // draw random year intercept

  for(n in 1:Nhold){
    mu_pred[n] <- exp(gint_out[n] + a_pred[yidhold[n] - nyrs ] + coverEff_pred[n]);
    lambda_pred[n] <- trueP1_pred[n, spp]*mu_pred[n];
  }

  for(n in 1:Nhold){
    log_lik2[n] <- neg_binomial_2_log(Yhold[n], lambda_pred[n], theta);
  }
  
}

