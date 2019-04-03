data{
  int<lower=0> N;                   // observations
  int<lower=0> Y[N];                // observation vector
  int<lower=0> K;                   // # of fixed effects
  matrix[N,K] X;                    // fixed effects matrix   
  int<lower=0> G;                   // groups
  int<lower=0> g[N];                // group id
  vector[N] P2intra; 
  
  // for out of sample prediction
  // int<lower=0> hold_N;              // observations
  // int<lower=0> hold_Y[hold_N];      // observation vector
  // 
  // matrix[hold_N,K] hold_X;          // fixed effects matrix   
  // row_vector[J] hold_Z[hold_N];     // group level effects matrix
  // int<lower=0> hold_G;              // groups
  // int<lower=0> hold_g[hold_N];      // group id
  // 
  // matrix[hold_N, Nspp] hold_parents1;   // parents in plot
  // matrix[hold_N, Nspp] hold_parents2;   // parents in group
  // 
  // // all data for IBM 
  // int<lower=0, upper=1> IBM;          // Flag if predictions for IBM are to be generated 
  // 
  // int<lower=0> IBM_N;              // observations
  // int<lower=0> IBM_Y[IBM_N];      // observation vector
  // 
  // matrix[IBM_N,K] IBM_X;          // fixed effects matrix   
  // row_vector[J] IBM_Z[IBM_N];     // group level effects matrix
  // int<lower=0> IBM_G;              // groups
  // int<lower=0> IBM_g[IBM_N];      // group id
  // 
  // matrix[IBM_N, Nspp] IBM_parents1;   // parents in plot
  // matrix[IBM_N, Nspp] IBM_parents2;   // parents in group
}
parameters{
  vector[K] beta;               // fixed effects 
	vector[G] u_raw;                  // raw group-level effects 
  real<lower=0> u_scale;
  real<lower=0> theta;          // negative binomial scale 
}
transformed parameters{
  vector[N] lambda;
  vector[N] fixef;
  vector[G] u; 
  
  fixef = X*beta;
  
  for( i in 1:G){
    u[i] = u_raw[i]*u_scale; 
  }
  
  for(n in 1:N){  
    lambda[n] = exp(log(P2intra[n]) + fixef[n] + u[g[n]]);
  }
}
model{
  // Priors
  beta ~ normal(0,1);
  u_raw ~ normal(0,1);
  theta ~ cauchy(0,1);
  u_scale ~ cauchy(0,1);

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
  
  for(n in 1:N){ 
    log_lik[n] = neg_binomial_2_lpmf( Y[n] | lambda[n], theta); 
    Y_hat[n] = neg_binomial_2_rng( lambda[n], theta);
  }
     
  // 1. Holdout data predictions 
  
  // if( hold_N > 1 ){ 
  //   
  //   vector[J] hold_u[hold_G];
  //   matrix[J, hold_G] hold_u_raw;
  //   
  //   for(i in 1:hold_G)
  //     for(j in 1:J)
  //       hold_u_raw[j, i] = normal_rng(0,1);
  // 
  //   for(j in 1:hold_G)
  //     hold_u[j] = Sigma_L * col(hold_u_raw, j);
  // 
  //   hold_fixef = hold_X*beta;
  //   hold_trueP1 = hold_parents1*m + hold_parents2*(1-m);
  //   
  //   for(n in 1:hold_N)
  //     for( j in 1:Nspp)
  //       hold_trueP2[n, j] = sqrt(hold_trueP1[n, j]);
  //   
  //   hold_coverEff = hold_trueP2*w;
  // 
  //   for(n in 1:hold_N){
  //     hold_mu[n] = exp(hold_fixef[n] + hold_coverEff[n] + hold_Z[n]*hold_u[hold_g[n]]);
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
  // 
  // // 2. all data for IBM  
  // if( IBM ==  1 ){ 
  //   vector[IBM_N] IBM_mu;       // linear predictor 
  //   vector[IBM_N] IBM_lambda;
  //   matrix[IBM_N, Nspp] IBM_trueP1;
  //   matrix[IBM_N, Nspp] IBM_trueP2;
  //   vector[J] IBM_u[IBM_G];
  //   matrix[J, IBM_G] IBM_u_raw;
  //   vector[IBM_N] IBM_fixef;
  //   vector[IBM_N] IBM_coverEff;
  // 
  //   for(i in 1:IBM_G)
  //     for(j in 1:J)
  //       IBM_u_raw[j, i] = normal_rng(0,1);
  // 
  //   for(j in 1:IBM_G)
  //     IBM_u[j] = Sigma_L * col(IBM_u_raw, j);
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
  //     IBM_mu[n] = exp(IBM_fixef[n] + IBM_coverEff[n] + IBM_Z[n]*IBM_u[IBM_g[n]]);
  //     IBM_lambda[n] = IBM_trueP1[n, spp]*IBM_mu[n];  
  //     IBM_Y_hat[n] = neg_binomial_2_rng( IBM_lambda[n], theta); 
  //   }
  // }else {
  //   IBM_Y_hat = to_vector(rep_array(0, IBM_N));
  // }
}

