data{
  int<lower=0> N;                   // observations
  int<lower=0> Y[N];                // observation vector
  
  int<lower=0> K;                   // # of fixed effects
  matrix[N,K] X;                    // fixed effects matrix with open space in plot as intercept 
  int<lower=0> G;                   // groups
  int<lower=0> g[N];                // group id
  
  vector[N] P1;                     // intraspecific cover in plot 
  vector[N] P2;                     // average intraspecific cover in plot group 

  // for out of sample prediction
  int<lower=0> hold_N;              // observations
  int<lower=0> hold_Y[hold_N];      // observation vector

  matrix[hold_N,K] hold_X;          // fixed effects matrix
  int<lower=0> hold_G;              // groups
  int<lower=0> hold_g[hold_N];      // group id

  vector[hold_N] hold_P1;                     // intraspecific cover in plot 
  vector[hold_N] hold_P2;                     // average intraspecific cover in plot group 

  // all data for IBM 
  int<lower=0, upper=1> IBM;          // Flag if predictions for IBM are to be generated

  int<lower=0> IBM_N;              // observations
  int<lower=0> IBM_Y[IBM_N];      // observation vector

  matrix[IBM_N,K] IBM_X;          // fixed effects matrix
  int<lower=0> IBM_G;              // groups
  int<lower=0> IBM_g[IBM_N];      // group id

  vector[IBM_N] IBM_P1;                     // intraspecific cover in plot 
  vector[IBM_N] IBM_P2;                     // average intraspecific cover in plot group 
}
parameters{
  vector[K] beta;                // fixed effects 
	real<lower=0> tau;             // scale parameter for group-level effects
	vector[G] u_raw;               // raw group-level effects
  real<lower=0> theta;           // negative binomial scale 
  real w;                        // intraspecific recruitment effect 
  real<lower=0, upper=1> m;      // mixing parameter for influence of outside cover
}
transformed parameters{
  vector[G] u;                // group-level effects
  vector[N] trueP1;
  vector[N] lambda;
  vector[N] fixef;

  for(j in 1:G)
    u[j] = u_raw[j] * tau;
  
  trueP1 = P1*m + P2*(1-m);
  
  fixef = X*beta;
    
  for(n in 1:N){
    lambda[n] = exp(fixef[n] + w*log(trueP1[n]) + u[g[n]]);
  }
}
model{
  // Priors
  beta ~ normal(0,1);
  tau ~ cauchy(0,1);                     // gamma as per rstanarm glmer vignette
  u_raw ~ normal(0,1);
  m ~ uniform(0,1);
  w ~ normal(0, 1);
  theta ~ cauchy(0, 1);

  // Likelihood
  Y ~ neg_binomial_2(lambda, theta);
}
generated quantities{
  vector[N] log_lik; 
  vector[N] Y_hat; 
  // vector[hold_N] hold_log_lik;
  // vector[hold_N] hold_mu;       // linear predictor
  // matrix[hold_N, Nspp] hold_trueP1;
  // matrix[hold_N, Nspp] hold_trueP2;
  // vector[hold_N] hold_lambda;
  // vector[hold_N] hold_coverEff;
  // vector[hold_N] hold_fixef;
  // vector[hold_N] hold_SE;
  // real hold_SSE;
  // vector[IBM_N] IBM_Y_hat;
  // 
  for(n in 1:N){ 
    log_lik[n] = neg_binomial_2_lpmf( Y[n] | lambda[n], theta); 
    Y_hat[n] = neg_binomial_2_rng( lambda[n], theta);
  }
     
  // 1. Holdout data predictions 
  
  // if( hold_N > 1 ){
  // 
  //   vector[hold_G] hold_u;
  // 
  //   for(j in 1:hold_G)
  //     hold_u[j] = normal_rng(0,1)*tau;
  // 
  //   hold_trueP1 = hold_parents1*m + hold_parents2*(1-m);
  // 
  //   for(n in 1:hold_N)
  //     for( j in 1:Nspp)
  //       hold_trueP2[n, j] = log(hold_trueP1[n, j]);
  // 
  //   hold_coverEff = hold_trueP2[,2]*w;
  //   
  //   hold_fixef = hold_X*beta;
  // 
  //   for(n in 1:hold_N){
  //     hold_mu[n] = exp(hold_fixef[n] + hold_coverEff[n] + hold_u[hold_g[n]]);
  //     hold_lambda[n] = hold_trueP1[n, spp]*hold_mu[n];
  //     hold_log_lik[n] = neg_binomial_2_lpmf( hold_Y[n] | hold_lambda[n], theta);
  //   }
  // }else if(hold_N <= 2 ){
  //     hold_fixef = to_vector(rep_array(0, hold_N));
  //     hold_mu = to_vector(rep_array(0, hold_N));
  //     hold_lambda = to_vector(rep_array(0, hold_N));
  //     hold_log_lik = to_vector(rep_array(negative_infinity(), hold_N));
  // }
  // 
  // for( i in 1:hold_N){
  //   hold_SE[i] = (hold_lambda[i] - hold_Y[i])^2 ;
  // }
  // hold_SSE = sum(hold_SE)/hold_N;
  // 

  // 2. all data for IBM  
  // if( IBM ==  1 ){
  //   vector[IBM_N] IBM_mu;       // linear predictor
  //   vector[IBM_N] IBM_lambda;
  //   matrix[IBM_N, Nspp] IBM_trueP1;
  //   matrix[IBM_N, Nspp] IBM_trueP2;
  //   vector[IBM_G] IBM_u;
  //   vector[IBM_N] IBM_fixef;
  //   vector[IBM_N] IBM_coverEff;
  // 
  //   for(j in 1:IBM_G)
  //     IBM_u[j] = normal_rng(0,1)*tau;
  // 
  //   IBM_fixef = IBM_X*beta;
  //   IBM_trueP1 = IBM_parents1*m + IBM_parents2*(1-m);
  // 
  //   for(n in 1:IBM_N)
  //     for( j in 1:Nspp)
  //       IBM_trueP2[n, j] = sqrt(IBM_trueP1[n, j]);
  // 
  //   IBM_coverEff = IBM_trueP2*w;
  // 
  //   for(n in 1:IBM_N){
  //     IBM_mu[n] = exp(IBM_fixef[n] + IBM_coverEff[n] + IBM_u[IBM_g[n]]);
  //     IBM_lambda[n] = IBM_trueP1[n, spp]*IBM_mu[n];
  //     IBM_Y_hat[n] = neg_binomial_2_rng( IBM_lambda[n], theta);
  //   }
  // }else {
  //   IBM_Y_hat = to_vector(rep_array(0, IBM_N));
  // }
}

