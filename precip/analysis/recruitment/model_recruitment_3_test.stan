// Single species model with climate: includes climate + intraspecific effects
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> Nspp; // number of species 
  int<lower=0> spp; // focal species id
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  
  int<lower=0> Covs; // climate covariates
  matrix[N,Covs] C; // climate matrix
  real tau_beta;
  
  
}parameters{
  real a_mu;
  vector[Yrs] a;
  vector[Nspp] w;
  vector[Covs] b2;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] q;
  vector[N] climEff;
  vector[N] coverEff;
  vector[N] p1; 
  vector[N] p2;
  vector[N] trueP1vec;
  vector[N] trueP2vec;
  
  p1 <- parents1[, Nspp];
  p2 <- parents2[, Nspp];

  trueP1 <- parents1*0.5 + parents2*(1-0.5);
  
  trueP1vec <- p1*0.5 + p2*(1-0.5); 
  
  climEff <- C*b2;

  for(n in 1:N){
    for( j in 1:Nspp) { 
      trueP1[n, j] <- parents1[n, j]*u + parents2[n,j]*(1-u);
      //print( trueP1[n,j])
      print( sqrt(trueP1[n, j]));
    }
    trueP2vec[n] <- sqrt(trueP1vec[n]);
  }
  
  coverEff <- trueP2vec*w[Nspp];

  for(n in 1:N){
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + coverEff[n] + climEff[n]);
    lambda[n] <- trueP1vec[n]*mu[n];  // elementwise multiplication  
  } 
  
  q <- lambda*theta;
  //print( trueP2vec[1:10])
  //print( trueP2[1:10, Nspp])
  //print(Nspp)
  
}
model{
  // Priors
  u ~ uniform(0.1,0.9);
  theta ~ cauchy(0,2);
  a_mu ~ normal(0,5);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  w ~ normal(0, 5);
  b2 ~ normal(0, tau_beta);
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
generated quantities{
  
  vector[N] log_lik; // vector for computing log pointwise predictive density
  
  for(n in 1:N)
    log_lik[n] <- neg_binomial_2_log(Y[n], q[n], theta); 

}

