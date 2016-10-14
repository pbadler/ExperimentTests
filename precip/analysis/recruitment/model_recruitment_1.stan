// Single-species model for recruitment: includes intraspecific effects 
data{
  int<lower=0> N; // observations  
  int<lower=0> Y[N]; // observation vector
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Nspp; // number of species 
  int<lower=0> spp; // focal species id
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  
}parameters{
  real a_mu;
  vector[Yrs] a;
  real w;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
  vector[N] mu;
  vector[N] trueP1;
  vector[N] trueP2;
  vector[N] lambda;
  vector[N] q;
  vector[N] coverEff;
  vector[N] p1; 
  vector[N] p2;
  
  p1 <- parents1[, spp];
  p2 <- parents2[, spp];

  trueP1 <- p1*u + p2*(1-u);

  for(n in 1:N)
      trueP2[n] <- sqrt(trueP1[n]);
  
  coverEff <- trueP2*w;

  for(n in 1:N){
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + coverEff[n]);
    lambda[n] <- trueP1[n]*mu[n];  // elementwise multiplication  
  } 
  
  q <- lambda*theta;

}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,5);
  a_mu ~ normal(0,5);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  w ~ normal(0, 5);
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
