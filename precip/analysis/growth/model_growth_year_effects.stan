data{
  // training datalist, historical observations 
  int<lower=0> N;              // observations
  vector[N] Y;                // observation vector
  vector[N] X;                // size vector
 int<lower=0> G;             // groups
  matrix[N, G] gm;            // group dummy variable matrix 
  #int<lower=0> gid[N];         // group id 
  int<lower=0> nyrs;           // years
  int<lower=0> yid[N];        // year id
  int<lower=0> Wcovs;           // number of crowding effects 
  matrix[N,Wcovs] W;          // crowding matrix
}
parameters{
  // for training data model  
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  vector[Wcovs] w;
  real<lower=0> sigma; 
  vector[G] bg;                     // varying group effects with first group as intercept 

}
transformed parameters{
  // for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff; 
  vector[N] gint; 

  // for training data model -----------------------------------
  crowdEff  <- W*w;
  gint <- gm*bg;

  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n]     <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n];
  }
}
model{
  // for training data model 
  // Priors
  sigma ~ cauchy(0,5);
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);
  
  // Likelihood
  Y ~ normal(mu, sigma);
}
