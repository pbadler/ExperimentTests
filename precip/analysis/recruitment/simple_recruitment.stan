// Full model for recruitment: includes climate + intraspecific + interspecific competition effects 
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  vector[N] parents1; // parents in plot
  vector[N] parents2; // parents in group
  
  int<lower=0> Covs; // climate covariates
  matrix[N,Covs] C; // climate matrix
  real tau_beta;
  
}parameters{
  real a_mu;
  vector[Yrs] a;
  real dd;
  vector[Covs] b2;

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
  vector[N] climEff;

  
  climEff <- C*b2;
  
  trueP1 <- parents1*u + parents2*(1-u);
  
  for( n in 1:N)
    trueP2[n] <- sqrt(trueP1[n]);
  
  for(n in 1:N)
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + dd*trueP2[n] + climEff[n]);
  
  lambda <- trueP1 .* mu;  // note use of elementwise multiplication operator 
    
  q <- lambda*theta;
    
}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,10);
  a_mu ~ normal(0,10);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  dd ~ uniform(-10, 10);
  
  b2 ~ normal(0, tau_beta);
  
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
