data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  vector[N] X;                // size vector
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  int<lower=0> Wcovs;         // number of crowding effects 
  vector[N] W;                // crowding matrix

  // survival data for cover predictions
  int<lower=0> N3;              // observations
  vector[N3] X3;                // size vector
  int<lower=0> yid3[N3];        // year id
  vector[N3] W3;                // crowding matrix
}
parameters{
  // for training data model  
  real<lower=0> sigma; 
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real w;
  real b0; 
}
transformed parameters{
// for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  real mu[N];
  vector[N] crowdEff;

  // for training data model -----------------------------------
  crowdEff <- W*w;

  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- b0 + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n];
  }
}
model{
   // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  b0 ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);

  // Likelihood
  Y ~ normal(mu, sigma);
  
}
generated quantities {
  // hold out predictions 
  real mu_out[N3];
  vector[N3] crowd_out;

  crowd_out   <- W3*w;

  for(n in 1:N3){
      mu_out[n] <- b0 + a[yid3[n]] + b1[yid3[n]]*X3[n] + crowd_out[n];
  }

}
