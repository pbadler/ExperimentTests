data{
  int<lower=0> N;             // observations
  int<lower=0> Y[N];          // observation vector
  
  int<lower=0> K;                       // # of fixed effects
  matrix[N,K] X;                        // fixed effects matrix   
  int<lower=0> J;                       // # of group level effects
  row_vector[J] Z[N];                   // group level effects matrix
  int<lower=0> G;                       // groups
  int<lower=0> g[N];                    // group id
  
  
  int<lower=0> Nspp;          // number of species 
  int<lower=0> spp;           // focal species id
  matrix[N, Nspp] parents1;   // parents in plot
  matrix[N, Nspp] parents2;   // parents in group

  // for out of sample prediction
  // int<lower=0> Nhold;
  // int<lower=0> nyrshold;              // years out
  // int<lower=0> yidhold[Nhold];        // year out id
  // int<lower=0> Yhold[Nhold];          // observation vector
  // matrix[Nhold, Nspp] parents1hold;   // hold out parents in plot
  // matrix[Nhold, Nspp] parents2hold;   // hold out parents in group
  // matrix[Nhold, G] gmhold;
  // matrix[Nhold,Covs] Chold;           // climate matrix
  
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
// generated quantities{
//   vector[N] log_lik; 
//   vector[Nhold] log_lik2; 
//   vector[nyrshold] a_pred;
//   vector[Nhold] mu_pred;
//   vector[Nhold] lambda_pred;
// 
//   // for prediction 
//   vector[Nhold] coverEff_pred;
//   matrix[Nhold, Nspp] trueP1_pred;
//   matrix[Nhold, Nspp] trueP2_pred;
//   vector[Nhold] gint_out; 
//   vector[Nhold] climhat; 
//   
//   // for cover prediction 
//   vector[N2] coverEff_pred2;
//   matrix[N2, Nspp] trueP1_pred2;
//   matrix[N2, Nspp] trueP2_pred2;
//   vector[N2] gint_out2; 
//   vector[N2] climhat2; 
//   vector[nyrs2] a_pred2;  
//   vector[N2] mu_pred2;
//   vector[N2] lambda_pred2;
//   
//   
//   for(n in 1:N){ 
//     log_lik[n] <- neg_binomial_2_log(Y[n], lambda[n], theta); 
//   }
//   
//   // 1. Holdout data predictions 
// 
//   climhat <- Chold*b2;
//   gint_out   <- gmhold*bg;
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
// }
