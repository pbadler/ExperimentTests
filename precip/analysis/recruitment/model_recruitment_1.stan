// Single-species model for recruitment: includes intraspecific effects 
data{
  int<lower=0> N; // observations  
  int<lower=0> Y[N]; // observation vector
  int<lower=0> nyrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Nspp; // number of species 
  int<lower=0> spp; // focal species id
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  matrix[N, G] gm;
}parameters{
  vector[nyrs] a_raw;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  real w;
  vector[G] bg; 
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
  vector[N] gint; 
  vector[nyrs] a; 
  
  for( n in 1:N){
    p1[n] <- parents1[n,spp];
    p2[n] <- parents2[n,spp];
  }

  trueP1 <- p1*u + p2*(1-u);

  for(n in 1:N)
      trueP2[n] <- sqrt(trueP1[n]);
  
  gint     <- gm*bg;
  coverEff <- trueP2*w;
  
  a  <- 0 + a_raw*sig_a; 
  
  for(n in 1:N){
    mu[n] <- exp(gint[n]  + a[yid[n]]  + coverEff[n]);
    lambda[n] <- trueP1[n]*mu[n];  
    q[n] <- lambda[n]*theta; 
  } 

}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,5);
  sig_a ~ cauchy(0,5);
  w ~ normal(0, 5);
  a_raw ~ normal(0, 1);
  bg ~ normal(0, 10);
  
  // Likelihood
  Y ~ neg_binomial_2(q, theta);
  
}
generated quantities{
  
  vector[N] log_lik; // vector for computing log pointwise predictive density
  
  for(n in 1:N)
    log_lik[n] <- neg_binomial_2_log(Y[n], q[n], theta); 
}
