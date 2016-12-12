data{
  int<lower=0> Nhold;                 // observations
  int<lower=0> Yhold[Nhold];          // observation vector
  int<lower=0> nyrshold;              // years
  int<lower=0> nyrs;
  int<lower=0> yidhold[Nhold];        // year id
  int<lower=0> G;                     // Groups 
  matrix[Nhold, G] gmhold;            // group dummy variable matrix
  int<lower=0> nT;               // groups
  matrix[Nhold, nT] tmhold;      // group dummy variable matrix
  int<lower=0> Nspp;                  // number of species 
  int<lower=0> spp;                   // focal species id
  matrix[Nhold, Nspp] parents1hold;   // parents in plot
  matrix[Nhold, Nspp] parents2hold;   // parents in group
  
  // all data for cover predictions 
  int<lower=0> N2;
  int<lower=0> nyrs2;              // years out
  int<lower=0> yid2[N2];        // year out id
  int<lower=0> Y2[N2];          // observation vector
  matrix[N2, Nspp] parents12;   // hold out parents in plot
  matrix[N2, Nspp] parents22;   // hold out parents in group
  matrix[N2, G] gm2;
  matrix[N2, nT] tm2;      // group dummy variable matrix
  
}parameters{
  vector[nyrshold] a_raw;
  vector[Nspp] w;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  vector[nT] bt; 
  vector[G] bg;
}
transformed parameters{
  vector[Nhold] mu;
  matrix[Nhold, Nspp] trueP1;
  matrix[Nhold, Nspp] trueP2;
  vector[Nhold] lambda;
  vector[Nhold] coverEff;
  vector[Nhold] treatEff;
  vector[nyrshold] a; 
  vector[Nhold] gint; 
  vector[Nhold] year_effect;
  
  trueP1 <- parents1hold*u + parents2hold*(1-u);
  
  for(n in 1:Nhold)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  gint     <- gmhold*bg;
  treatEff <- tmhold*bt;
  coverEff <- trueP2*w;
  a  <- 0 + a_raw*sig_a; 

  for(n in 1:Nhold){
    year_effect[n] <- a[yidhold[n] - nyrs]; 
    
    mu[n] <- exp(treatEff[n] + year_effect[n] + coverEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  
  }
}

model{
  // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,2);
  sig_a ~ cauchy(0,2);
  a_raw ~ normal(0, 1);
  bt ~ normal(0, 10);
  bg ~ normal(0, 10);
  w ~ normal(0, 10);
  
  // Likelihood
  Yhold ~ neg_binomial_2(lambda, theta);
}
generated quantities{
  vector[Nhold] log_lik2; 
  
  // for cover prediction 
  vector[N2] coverEff_pred2;
  matrix[N2, Nspp] trueP1_pred2;
  matrix[N2, Nspp] trueP2_pred2;
  vector[N2] gint_out2; 
  vector[N2] treatEff2; 
  vector[nyrs2] a_pred2;  
  vector[N2] mu_pred2;
  vector[N2] lambda_pred2;
  
  for(n in 1:Nhold){ 
    log_lik2[n] <- neg_binomial_2_log(Yhold[n], lambda[n], theta); 
  }
  
  // 2. all data for cover predictions 
  treatEff2      <- tm2*bt;
  gint_out2     <- gm2*bg;
  trueP1_pred2  <- parents12*u + parents22*(1-u);

  for(n in 1:N2)
    for(j in 1:Nspp)
      trueP2_pred2[n, j] <- sqrt(trueP1_pred2[n, j]);
  
  coverEff_pred2 <- trueP2_pred2*w;

  for( i in 1:nyrs2)
    a_pred2[i] <- normal_rng(0, sig_a); // draw random year intercept

  for(n in 1:N2){
    mu_pred2[n] <- exp(gint_out2[n] + a_pred2[yid2[n]] + coverEff_pred2[n] + treatEff2[n]);
    lambda_pred2[n] <- trueP1_pred2[n, spp]*mu_pred2[n];
  }

}
