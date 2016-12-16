data{
  int<lower=0> N2;                 // observations
  int<lower=0> Y2[N2];          // observation vector
  int<lower=0> nyrs2;              // years
  int<lower=0> yid2[N2];        // year id
  int<lower=0> G;                     // Groups 
  matrix[N2, G] gm2;            // group dummy variable matrix
  int<lower=0> nT;               // groups
  matrix[N2, nT] tm2;      // group dummy variable matrix
  int<lower=0> Nspp;                  // number of species 
  int<lower=0> spp;                   // focal species id
  matrix[N2, Nspp] parents12;   // parents in plot
  matrix[N2, Nspp] parents22;   // parents in group
  
}parameters{
  vector[nyrs2] a_raw;
  vector[Nspp] w;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  vector[nT] bt; 
  vector[G] bg;
}
transformed parameters{
  vector[N2] mu;
  matrix[N2, Nspp] trueP1;
  matrix[N2, Nspp] trueP2;
  vector[N2] lambda;
  vector[N2] coverEff;
  vector[N2] treatEff;
  vector[nyrs2] a; 
  vector[N2] gint; 
  vector[N2] year_effect;
  
  trueP1 <- parents12*u + parents22*(1-u);
  
  for(n in 1:N2)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  gint     <- gm2*bg;
  treatEff <- tm2*bt;
  coverEff <- trueP2*w;
  a  <- 0 + a_raw*sig_a; 

  for(n in 1:N2){
    year_effect[n] <- a[yid2[n]]; 
    
    mu[n] <- exp(gint[n] + year_effect[n] + coverEff[n] + treatEff[n]);
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
  Y2 ~ neg_binomial_2(lambda, theta);
}
generated quantities{
  vector[N2] log_lik2; 
  
  // for cover prediction 
  vector[nyrs2] a_pred2;  
  vector[N2] mu_pred2;
  vector[N2] lambda_pred2;
  
  for(n in 1:N2){ 
    log_lik2[n] <- neg_binomial_2_log(Y2[n], lambda[n], theta); 
  }
  
  // 2. all data for cover predictions 
  for( i in 1:nyrs2)
    a_pred2[i] <- normal_rng(0, sig_a); // draw random year intercept

  for(n in 1:N2){
    mu_pred2[n] <- exp(gint[n] + a_pred2[yid2[n]] + coverEff[n] + treatEff[n]);
    lambda_pred2[n] <- trueP1[n, spp]*mu_pred2[n];
  }

}
