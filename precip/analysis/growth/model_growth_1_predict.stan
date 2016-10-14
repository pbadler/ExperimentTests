// Single species growth model: includes intraspecific effects only 
data{
  int<lower=0> N; // observations
  vector[N] Y; // observation vector
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  vector[N] X; // size vector
  int<lower=0> Wcovs; // number of crowding effects 
  matrix[N,Wcovs] W; // crowding matrix
  int<lower=0>Covs; // number of climate effects 
  matrix[N,Covs] C; // climate matrix
  real tau_beta;

  // For out of sample prediction
  int<lower=0> npreds;
  vector[npreds] y_holdout;
  vector[npreds] Xhold;
  int<lower=0> nyrs_out; // years holdout 
  int<lower=0> yid_out[npreds];// year id holdout
  int<lower=0> gid_out[npreds]; // group id holdout
  matrix[npreds,Wcovs] Whold; // crowding matrix for holdout data 
  matrix[npreds,Covs] Chold;
  
  // for year effect estimation using entire dataset 
  int<lower=0> N2; // all observations
  vector[N2] Y2; // observation vector
  int<lower=0> Yrs2; // all years
  int<lower=0> yid2[N2]; // year id
  int<lower=0> gid2[N2]; // group id 
  vector[N2] X2; // size vector
  matrix[N2,Wcovs] W2; // crowding matrix

  // For predicting overall cover 
  int<lower=0> N3; // observations
  int<lower=0> nyrs3; // years
  int<lower=0> yid3[N3]; // year id
  int<lower=0> gid3[N3]; // group id
  vector[N3] X3; // size vector
  matrix[N3,Wcovs] W3; // crowding matrix
  matrix[N3,Covs] C3; // climate matrix
  
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Wcovs] w;
  real gint[G];
  real tau;
  real tauSize;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
  
  // for year effects model  
  real a_mu2;
  vector[Yrs2] a2;
  real b1_mu2;
  vector[Yrs2] b12;
  real<lower=0> sig_a2;
  real<lower=0> sig_b12;
  real gint2[G];
  real tau2;
  real tauSize2;
  vector[Wcovs] w2;
  real<lower=0> sig_G2;
}
transformed parameters{
  vector[N] mu;
  real<lower=0> sigma[N];
  vector[N] crowdEff;
  vector[N] climEff;
  real mu2[N2] ;  
  real<lower=0> sigma2[N2];
  vector[N2] crowdEff2; 
  
  crowdEff <- W*w;

  for(n in 1:N){
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
  }
  
  // for year effects model 
  
  crowdEff2 <- W2*w2;
  
  for(n in 1:N2){
    mu2[n] <- a2[yid2[n]] + gint2[gid2[n]] + b12[yid2[n]]*X2[n] + crowdEff2[n];
    sigma2[n] <- sqrt((fmax(tau2*exp(tauSize2*mu2[n]), 0.0000001)));  
  }
  
}
model{
  // Priors
  a_mu ~ normal(0,10);
  b1_mu ~ normal(0,10);
  w ~ normal(0, 10);
  tau ~ normal(0,10);
  tauSize ~ normal(0,10);
  sig_a ~ cauchy(0,2);
  sig_b1 ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  for(g in 1:G)
    gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, sigma);
  
  //for year effects model 
  a_mu2 ~ normal(0,10);
  w2 ~ normal(0,10);
  b1_mu2 ~ normal(0,10);
  sig_a2 ~ cauchy(0,2);
  sig_b12 ~ cauchy(0,2);
  sig_G2 ~ cauchy(0,2);
  tau2 ~ normal(0,10);
  tauSize2 ~ normal(0,10);
  for(g in 1:G)
    gint2[g] ~ normal(0, sig_G2);
  for(y in 1:Yrs){
    a2[y] ~ normal(a_mu2, sig_a2);
    b12[y] ~ normal(b1_mu2, sig_b12);
  }

  // Likelihood
  Y2 ~ normal(mu2,sigma2);

}
generated quantities {
  vector[nyrs_out] a_out;
  vector[nyrs_out] b1_out;
  vector[npreds] crowdhat;
  vector[npreds] sigmahat;
  vector[npreds] muhat;
  vector[npreds] y_hat;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  
  // for predictions from known year effects model 
  vector[npreds] muhat2;
  vector[npreds] y_hat2; // pointwise predictions  
  vector[npreds] log_lik2; // vector for computing log pointwise predictive density  
  int<lower=0> yid_out2[npreds]; //integer for modern year effects  
  vector[npreds] sigmahat2; // standard deviation for modern predictions 
  
  // for cover predictions 
  vector[nyrs3] a3;
  vector[nyrs3] b13;
  vector[N3] crowdhat3;
  vector[N3] muhat3;
  vector[N3] y_hat3; // pointwise predictions  
  vector[N3] sigmahat3; // standard deviation for modern predictions 
  
  // for cover predictions from known year effects model
  int<lower=0>  yid4[N3]; // integer for year effects 
  vector[N3] muhat4;
  vector[N3] sigmahat4;
  vector[N3] y_hat4;
  
  // 1. Hold out data predictions 
  crowdhat <- Whold*w;
  
  for( i in 1:nyrs_out){
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1); //draw random year x size effect 
  }
  
  for(n in 1:npreds){
    muhat[n] <- a_out[yid_out[n]] + gint[gid_out[n]] + b1_out[yid_out[n]]*Xhold[n] + crowdhat[n];
    sigmahat[n] <- sqrt((fmax(tau*exp(tauSize*muhat[n]), 0.0000001))); 
    y_hat[n] <- normal_rng(muhat[n], sigmahat[n]);
    log_lik[n] <- normal_log(y_holdout[n], muhat[n], sigmahat[n]);
  }
  
  // 2. Predictions for holdout data with KNOWN year effects.  
  //    Simulate predictions as if year effects in the out of sample data are known. 
  
  for( n in 1:npreds){ 
    yid_out2[n] <- yid_out[n] + Yrs;  // add number of training years to get correct index for a2 and b12
    muhat2[n] <- a2[yid_out2[n]] + gint[gid_out[n]] + b12[yid_out2[n]]*Xhold[n] + crowdhat[n];
    sigmahat2[n] <- sqrt((fmax(tau*exp(tauSize*muhat2[n]), 0.0000001))); 
    y_hat2[n] <- normal_rng(muhat2[n], sigmahat2[n]);
    log_lik2[n] <- normal_log(y_holdout[n], muhat2[n], sigmahat2[n]);
  }
  
  // 3. Generate predicted size for all plants in the survival dataset 
  
  crowdhat3 <- W3*w;
  
  for( i in 1:nyrs3){
    a3[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 
    b13[i] <- normal_rng(b1_mu, sig_b1); //draw random year x size effect 
  }

  for( n in 1:N3){
    muhat3[n] <- a3[yid3[n]] + gint[gid3[n]] + b13[yid3[n]]*X3[n] + crowdhat3[n];
    sigmahat3[n] <- sqrt((fmax(tau*exp(tauSize*muhat3[n]), 0.0000001))); 
    y_hat3[n] <- normal_rng(muhat3[n], sigmahat3[n]);
  }
  
  // 4. Generate predicted size for all plants in the survival dataset 
  //    using KNOWN year effects. 
  
  for( n in 1:N3){
      yid4[n] <- yid3[n] + Yrs;  // add number of training years to get correct index for a4 and b14
      muhat4[n] <- a2[yid4[n]] + gint[gid3[n]] + b12[yid4[n]]*X3[n] + crowdhat3[n];
      sigmahat4[n] <- sqrt((fmax(tau*exp(tauSize*muhat4[n]), 0.0000001))); 
      y_hat4[n] <- normal_rng(muhat4[n], sigmahat4[n]);
  }
      
}

