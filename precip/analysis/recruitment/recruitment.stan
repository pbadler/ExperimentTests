data{
  int<lower=0> N;                   // observations
  int<lower=0> Y[N];                // observation vector
  
  int<lower=0> K;                   // # of fixed effects
  matrix[N,K] X;                    // fixed effects matrix   
  int<lower=0> J;                   // # of group level effects
  row_vector[J] Z[N];               // group level effects matrix
  int<lower=0> G;                   // groups
  int<lower=0> g[N];                // group id
  
  int<lower=0> Nspp;                // number of species 
  int<lower=0> spp;                 // focal species id
  matrix[N, Nspp] parents1;         // parents in plot
  matrix[N, Nspp] parents2;         // parents in group

  // for out of sample prediction
  int<lower=0> hold_N;              // observations
  int<lower=0> hold_Y[hold_N];      // observation vector
  
  matrix[hold_N,K] hold_X;          // fixed effects matrix   
  row_vector[J] hold_Z[hold_N];     // group level effects matrix
  int<lower=0> hold_G;              // groups
  int<lower=0> hold_g[hold_N];      // group id
  
  matrix[hold_N, Nspp] hold_parents1;   // parents in plot
  matrix[hold_N, Nspp] hold_parents2;   // parents in group

  // all data for cover predictions 
  // int<lower=0> N2;
  // int<lower=0> nyrs2;              // years out
  // int<lower=0> yid2[N2];        // year out id
  // int<lower=0> Y2[N2];          // observation vector
  // matrix[N2, Nspp] parents12;   // hold out parents in plot
  // matrix[N2, Nspp] parents22;   // hold out parents in group
  // matrix[N2, G] gm2;
  // matrix[N2,Covs] C2;           // climate matrix
  
}
transformed data{ 
  vector[J] a;                          // prior on dirichlet distribution

  for(i in 1:J)
    a[i] = 1; 
}
parameters{
  vector[K] beta;               // fixed effects 
  simplex[J] pi_;                // simplex for diagonal of group-level covariance matrix 
	real<lower=0> tau;             // scale parameter for group-level covariance matrix
  cholesky_factor_corr[J] L_u;   // cholesky factor for group-level correlation
	matrix[J,G] u_raw;             // raw group-level effects 

  real<lower=0> theta;           // negative binomial scale 
  vector[Nspp] w;
  real<lower=0, upper=1> m;      // mixing parameter for influence of outside cover
}
transformed parameters{
  vector[N] mu;               // linear predictor 
  vector[J] u[G];             // group-level effects 
  matrix[J,J] Sigma_L;        // cholesky of covariance matrix
  vector[J] sigma_j;          // diagonal of covariance matrix 
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] coverEff;
  vector[N] fixef;

  sigma_j = pi_*J*tau^2;                      // Explained by rstanarm glmer vignette
  Sigma_L = diag_pre_multiply(sigma_j, L_u);  // multiply variance by correlation  
  
  for(j in 1:G)
    u[j] = Sigma_L * col(u_raw, j);    

  fixef = X*beta;
  trueP1 = parents1*m + parents2*(1-m);
  
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] = sqrt(trueP1[n, j]);
  
  coverEff = trueP2*w;
    
  for(n in 1:N){
    mu[n] = exp(fixef[n] + coverEff[n] + Z[n]*u[g[n]]);
    lambda[n] = trueP1[n, spp]*mu[n];  
  }
}
model{
  // Priors
  beta ~ normal(0,5);
  pi_ ~ dirichlet(a);                   // dirichlet as per rstanarm glmer vignette
  tau ~ gamma(1,1);                     // gamma as per rstanarm glmer vignette
  L_u ~ lkj_corr_cholesky(1.0);
  to_vector(u_raw) ~ normal(0,1);
  m ~ uniform(0,1);
  w ~ normal(0, 5);
  theta ~ cauchy(0,2);

  // Likelihood
  Y ~ neg_binomial_2(lambda, theta);
}
generated quantities{
  vector[N] log_lik; 
  vector[hold_N] hold_log_lik; 
  vector[hold_N] hold_mu;       // linear predictor 
  matrix[hold_N, Nspp] hold_trueP1;
  matrix[hold_N, Nspp] hold_trueP2;
  vector[hold_N] hold_lambda;
  vector[hold_N] hold_coverEff;
  vector[hold_N] hold_fixef;
  
  vector[hold_N] hold_SE; 
  real hold_SSE; 
  
  for(n in 1:N){ 
    log_lik[n] = neg_binomial_2_lpmf( Y[n] | lambda[n], theta); 
  }
     
  // 1. Holdout data predictions 
  
  if( hold_N > 1 ){ 
    
    vector[J] hold_u[hold_G];
    matrix[J, hold_G] hold_u_raw;
    
    for(i in 1:hold_G)
      for(j in 1:J)
        hold_u_raw[j, i] = normal_rng(0,1);
  
    for(j in 1:hold_G)
      hold_u[j] = Sigma_L * col(hold_u_raw, j);
  
    hold_fixef = hold_X*beta;
    hold_trueP1 = hold_parents1*m + hold_parents2*(1-m);
    
    for(n in 1:hold_N)
      for( j in 1:Nspp)
        hold_trueP2[n, j] = sqrt(hold_trueP1[n, j]);
    
    hold_coverEff = hold_trueP2*w;
  
    for(n in 1:hold_N){
      hold_mu[n] = exp(hold_fixef[n] + hold_coverEff[n] + hold_Z[n]*hold_u[hold_g[n]]);
      hold_lambda[n] = hold_trueP1[n, spp]*hold_mu[n];  
      hold_log_lik[n] = neg_binomial_2_lpmf( hold_Y[n] | hold_lambda[n], theta); 
    }
  }else if(hold_N <= 2 ){
      hold_fixef = to_vector(rep_array(0, hold_N));
      hold_mu = to_vector(rep_array(0, hold_N));
      hold_lambda = to_vector(rep_array(0, hold_N));
      hold_log_lik = to_vector(rep_array(negative_infinity(), hold_N));
  }
  
  for( i in 1:hold_N){ 
    hold_SE[i] = (hold_lambda[i] - hold_Y[i])^2 ;
  }
  hold_SSE = sum(hold_SE)/hold_N;

//   trueP1_pred <- parents1hold*u + parents2hold*(1-u);
// 
//   for(n in 1:Nhold)
//     for(j in 1:Nspp)
//       trueP2_pred[n, j] <- sqrt(trueP1_pred[n, j]);
//   
//   coverEff_pred <- trueP2_pred*w;
// 
//   for( i in 1:nyrshold)
//     a_pred[i] <- normal_rng(0, sig_a); // draw random year intercept
// 
//   for(n in 1:Nhold){
//     mu_pred[n] <- exp(gint_out[n] + a_pred[yidhold[n] - nyrs ] + coverEff_pred[n] + climhat[n]);
//     lambda_pred[n] <- trueP1_pred[n, spp]*mu_pred[n];
//   }
// 
//   for(n in 1:Nhold){
//     log_lik2[n] <- neg_binomial_2_log(Yhold[n], lambda_pred[n], theta);
//   }
//   
//   // 2. all data for cover predictions 
//   climhat2      <- C2*b2;
//   gint_out2     <- gm2*bg;
//   trueP1_pred2  <- parents12*u + parents22*(1-u);
// 
//   for(n in 1:N2)
//     for(j in 1:Nspp)
//       trueP2_pred2[n, j] <- sqrt(trueP1_pred2[n, j]);
//   
//   coverEff_pred2 <- trueP2_pred2*w;
// 
//   for( i in 1:nyrs2)
//     a_pred2[i] <- normal_rng(0, sig_a); // draw random year intercept
// 
//   for(n in 1:N2){
//     mu_pred2[n] <- exp(gint_out2[n] + a_pred2[yid2[n]] + coverEff_pred2[n] + climhat2[n]);
//     lambda_pred2[n] <- trueP1_pred2[n, spp]*mu_pred2[n];
//   }
// 
//   
}

