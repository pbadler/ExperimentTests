data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  vector[N] X;                // size vector
  int<lower=0> G;             // groups
  matrix[N, G] gm;            // group dummy variable matrix 
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  int<lower=0> Wcovs;         // number of crowding effects 
  vector[N] W;                // crowding matrix
  int<lower=0> spp;           // focal species 
  int<lower=0> Covs;          // climate effects 
  matrix[N, Covs] C;          // climate matrix 
  real <lower=0> tau_beta;    // prior sd for climate effects 
  
  // full datalist, all observations for year effects  
  int<lower=0> N2;              // observations
  vector[N2] Y2;                // observation vector
  vector[N2] X2;                // size vector
  matrix[N2, G] gm2;            // group dummy variable matrix 
  int<lower=0> nyrs2;           // years
  int<lower=0> yid2[N2];        // year id
  vector[N2] W2;                // crowding matrix

  // holdout datalist, modern observations 
  int<lower=0> Nhold;
  vector[Nhold] Yhold;
  int<lower=0> nyrshold;            // years out
  int<lower=0> yidhold[Nhold];      //year out id
  matrix[Nhold, G] gmhold;          // group dummy variable matrix 
  vector[Nhold] Xhold;
  vector[Nhold] Whold;              // crowding matrix for holdout data
  matrix[Nhold, Covs] Chold;        // climate matrix 
  
  // survival data for cover predictions
  int<lower=0> N3;              // observations
  vector[N3] Y3;                // observation vector
  vector[N3] X3;                // size vector
  matrix[N3, G] gm3;            // group dummy variable matrix 
  int<lower=0> nyrs3;           // years
  int<lower=0> yid3[N3];        // year id
  vector[N3] W3;                // crowding matrix
  matrix[N3, Covs] C3;          // climate matrix 
}
parameters{
  // for training data model  
  real<lower=0> sigma; 
  vector[G] bg;                     // varying group effects with first group as intercept 
  vector[nyrs] a_raw;
  real b1_mu;
  vector[nyrs] b1_raw;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real w;
  vector[Covs]  b2; 
  
  // for year effects model 
  real<lower=0> sigma2;
  vector[G] bg2;                     // varying group effects with first group as intercept 
  vector[nyrs2] a_raw2;
  real b1_mu2;
  vector[nyrs2] b1_raw2;
  real<lower=0> sig_a2;
  real<lower=0> sig_b12;
  real w2;
}
transformed parameters{
// for training data model  
  vector[nyrs] a;
  vector[nyrs] b1;
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff; 
  
  // for year effects model 
  vector[nyrs2] a2;
  vector[nyrs2] b12;
  vector[N2] gint2; 
  real mu2[N2];
  vector[N2] crowdEff2;

  // for training data model -----------------------------------
  gint <- gm*bg;
  crowdEff <- W*w;
  climEff  <- C*b2;
  
  b1 <- b1_mu + sig_b1*b1_raw;
  a  <- 0 + sig_a*a_raw; 
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n];
  }
  
  // for year effects model -----------------------------------
  gint2 <- gm2*bg2;
  crowdEff2 <- W2*w2;
  
  b12 <- b1_mu2 + sig_b12*b1_raw2;
  a2  <- 0 + sig_a2*a_raw2; 
  
  for(n in 1:N2){
    mu2[n] <- gint2[n] + a2[yid2[n]] + b12[yid2[n]]*X2[n] + crowdEff2[n];
  }
  
}
model{
   // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  w ~ normal(0,10);
  b2 ~ normal(0,tau_beta); 
    
  // Likelihood
  Y ~ normal(mu, sigma);
  
  // for year effects model 
  // Priors
  sigma2 ~ cauchy(0, 5);
  bg2 ~ normal(0,10);
  b1_mu2 ~ normal(0,10);
  sig_a2 ~ cauchy(0,5);
  sig_b12 ~ cauchy(0,5);
  a_raw2 ~ normal(0,1);
  b1_raw2 ~ normal(0,1);
  w2 ~ normal(0,10);

  // Likelihood
  Y2 ~ normal(mu2, sigma2);
}
generated quantities {
  // hold out predictions 
  vector[nyrshold] a_out;
  vector[nyrshold] b1_out;
  real muhat[Nhold];
  vector[Nhold] gint_out;
  vector[Nhold] y_hat;                        // pointwise predictions  
  vector[Nhold] log_lik;                      // vector for computing log pointwise predictive density  
  vector[Nhold] crowdhat;
  vector[Nhold] climhat;
  
  // for predictions with year effects 
  real muhat2[Nhold];
  vector[Nhold] y_hat2;                        // pointwise predictions  
  vector[Nhold] log_lik2;                      // vector for computing log pointwise predictive density  
  
  // for cover predictions 
  vector[N3] crowdhat3;
  vector[N3] climhat3;
  vector[N3] gint3;
  real muhat3[N3];
  vector[N3] y_hat3;                          // pointwise predictions  
  
  // for cover predictions 
  real muhat4[N3];
  vector[N3] y_hat4;                          // pointwise predictions  
  
  // 1. Holdout data predictions 
  gint_out   <- gmhold*bg; 
  crowdhat <- Whold*w;
  climhat  <- Chold*b2; 
  
  for( i in 1:nyrshold){
    a_out[i] <- normal_rng(0, sig_a);         // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1);   //draw random year x size effect 
  }
  
  for(n in 1:Nhold){
      muhat[n] <- gint_out[n] + a_out[yidhold[n]-nyrs] + b1_out[yidhold[n]-nyrs]*Xhold[n] + crowdhat[n] + climhat[n];
      y_hat[n] <- normal_rng(muhat[n], sigma);
      log_lik[n] <- normal_log(Yhold[n], muhat[n], sigma);
  }

  // 2. Holdout data predictions with known year effects 
  
  for(n in 1:Nhold){
      muhat2[n] <- gint_out[n] + a2[yidhold[n]] + b12[yidhold[n]]*Xhold[n] + crowdhat[n];
      y_hat2[n] <- normal_rng(muhat2[n], sigma);
      log_lik2[n] <- normal_log(Yhold[n], muhat2[n], sigma);
  }
  
  // 3. Holdout survival data for predicting cover 
  gint3     <- gm3*bg; 
  crowdhat3 <- W3*w;
  climhat3  <- C3*b2; 
  
  for(n in 1:N3){
      muhat3[n] <- gint3[n] + a_out[yid3[n]-nyrs] + b1_out[yid3[n]-nyrs]*X3[n] + crowdhat3[n] + climhat3[n];
      y_hat3[n] <- normal_rng(muhat3[n], sigma);
  }
  
  // 4. Holdout survival data for predicting cover 
  for(n in 1:N3){
      muhat4[n] <- gint3[n] + a2[yid3[n]] + b12[yid3[n]]*X3[n] + crowdhat3[n];
      y_hat4[n] <- normal_rng(muhat3[n], sigma);
  }
}
