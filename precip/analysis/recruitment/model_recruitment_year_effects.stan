data{
  int<lower=0> N;                 // observations
  int<lower=0> Y[N];              // observation vector
  int<lower=0> nyrs;              // years
  int<lower=0> yid[N];            // year id
  int<lower=0> G;                 // Groups 
  matrix[N, G] gm;                // group dummy variable matrix
  int<lower=0> Nspp;              // number of species 
  int<lower=0> spp;               // focal species id
  matrix[N, Nspp] parents1;       // parents in plot
  matrix[N, Nspp] parents2;       // parents in group
  
}parameters{
  vector[nyrs] a_raw;
  vector[Nspp] w;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  vector[G] bg;
}
transformed parameters{
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] coverEff;
  vector[N] treatEff;
  vector[nyrs] a; 
  vector[N] gint; 
  
  trueP1 <- parents1*u + parents2*(1-u);
  
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  gint     <- gm*bg;
  coverEff <- trueP2*w;
  a  <- 0 + a_raw*sig_a; 

  for(n in 1:N){
    mu[n] <- exp(gint[n] + a[yid[n]] + coverEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  
  }
}

model{
  // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,2);
  sig_a ~ cauchy(0,2);
  a_raw ~ normal(0, 1);
  bg ~ normal(0, 10);
  w ~ normal(0, 10);
  
  // Likelihood
  Y ~ neg_binomial_2(lambda, theta);
}

