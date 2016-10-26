// Single-species model for recruitment: includes intraspecific effects 
data{
  int<lower=0> N; // observations  
  int<lower=0> Y[N]; // observation vector
  vector[N] X; // covariate 
}parameters{
  real beta0; 
  real beta1; 
  real<lower=0> phi;
}
transformed parameters{
  vector[N] mu;

  mu <- exp(beta0 + beta1*X); 

}
model{
  // Priors
  phi ~ cauchy(0,5);
  beta0 ~ normal(0, 10); 
  beta1 ~ normal(0, 10); 
  
  // Likelihood
  Y ~ neg_binomial_2(mu, phi);
}

