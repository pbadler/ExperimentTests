extract_data <- function(df){ 
  N <- nrow(df)
  X <- df$X
  K <- ncol(X)
  Z <- df$Z 
  J <- ncol(Z)
  g <- as.numeric(factor(df$g))
  G <- length(unique(g))
  Y <- as.numeric(df$Y)
  S <- as.numeric(df$survives)
  E <- df$E
  D <- ncol(E)
  
  obs <- which(df$censored == 0)
  cens <- which(df$censored == 1)
  censored <- as.numeric(df$censored)
  N_obs <- length(obs)
  N_cens <- length(cens)
  Y_obs  <- Y[obs]
  U <- unique( df$U )
  
  parents1 <- df$parents1
  parents2 <- df$parents2
  Nspp <- ncol ( df$parents1 )
  spp <- unique( df$spp ) 
  years <- unique( df$year )
  
  rm(df)
  out <- lapply( ls(), function(x) eval(parse(text = x)))
  names(out) <- ls()[!ls() == 'out']
  
  return(out)
}


split_df <- function(df, hold){ 
  if(all(hold == 0)){
    df_out <- split(df, df$g %in% hold)
    df_out$True <- data.frame(Y = rep(0,2))
    df_out$True$survives = rep( 0, 2)
    df_out$True$X = matrix(0, ncol = ncol(df$X), nrow = 2)
    df_out$True$Z = matrix(0, ncol = ncol(df$Z), nrow = 2)
    df_out$True$g = rep(0,2)
    df_out$True$E = matrix(0, ncol = ncol(df$E), nrow = 2)

  }else if(any(hold > 0)){ 
    df_out <- split(df, df$g %in% hold)
  }
  names(df_out) <- c('train', 'hold')
  return(df_out)
}

make_dl <- function(df){ 
  dl <- unlist( lapply( df, extract_data), recursive = F)
  names(dl) <- gsub( 'train.', '', names(dl) )
  names(dl) <- gsub( 'hold.', 'hold_', names(dl))
  return(dl)
}

get_lpd <- function(my_fit){ 
  require(loo)
  
  lpd <- log(colMeans(exp(extract_log_lik(my_fit, 'hold_log_lik'))))

  return(lpd)
}

plot_x_y <- function(myfit, X, Y, iter, bt = F){
  
  if(bt){
    Y_hat <- extract( myfit, 'Y_hat_bt')$Y_hat_bt
  }else if(!bt){
    Y_hat <- extract( myfit, 'Y_hat')$Y_hat
  }
  
  ylims = c(min(c(Y_hat[iter,], Y)), max(c(Y_hat[iter, ], Y)))
  
  par(mfrow = c(2,2))
  
  plot(X, Y, ylim = ylims)
  abline(0,1)
  
  plot(X, Y_hat[iter,], ylim = ylims)
  abline(0,1)
  
  plot(Y, Y_hat[iter, ], ylim = ylims)
  abline(0,1)
  
  par(mfrow = c(1,1))
}

left_censor_df <- function(df, left_cut){ 
  # account for left censored data 
  # left cut is the cut off for the censored data on the original scale
  # of logarea.t1
  # it is rescaled to the Y-scale and then applied 
  
  U <- scale(left_cut,  attributes(df$Y)$`scaled:center`, attributes(df$Y)$`scaled:scale` )
  print(U)
  df$U <- as.numeric(U)
  df$censored <- as.numeric( df$Y <= df$U )
 
  return(df)
}


init_norm <- function(means=0, sds=1, name = NULL){
  out <- list(mapply( x = means, y = sds, function(x, y) rnorm(1, x, y)))
  names(out) <- name 
  out 
}

init_unif <- function(lower, upper, name = NULL){
  out <- list(mapply( x = lower, y = upper, function(x, y) runif(1, x, y)))
  names(out) <- name 
  out 
}


inv_logit <- function(x) { exp(x)/(1 + exp(x)) }

init_chol <- function(r_mu, r_sd, nr = 2, nc = 2, name = NULL){
  r <- 2*inv_logit(rnorm(1, r_mu, r_sd)) - 1
  out <- list(t(chol(matrix(c(1,r,r,1), 2,2))))
  names(out) <- name 
  out 
}

get_spp_and_vr <- function(dat_file, model_file){ 
  spp <- unlist( str_extract_all(c('ARTR', 'HECO', 'POSE' , 'PSSP'), string = dat_file) )
  vr <- unlist(str_extract_all(c('growth', 'survival', 'recruitment'), string = model_file))
  return(list(spp, vr))
}

process_data <- function(dat, formX, formC, formZ, formE = as.formula(~ -1), vr = 'growth', ... ){
  
  C <- model.matrix(formC, dat)
  dat$C <- scale(C)
  dat$W <- scale(dat$W)
  dat$Group <- factor(dat$gid)
  
  if( ncol(dat$C) == 0 ){ 
    formX <- update(formX, ' ~ . - C')
  }
  
  dat$X <- model.matrix(formX, data = dat)
  dat$Z <- model.matrix(formZ, data = dat)
  dat$E <- model.matrix(formE, data = dat) 
  
  dat$g <- factor(dat$yid)
  
  dat_4_cover <- dat ### Need to preserve dataframe with NA's (dead plants) for predicting cover 
  dat_4_cover <- split_df(dat_4_cover, hold = 0)
  dl_4_cover <- make_dl(dat_4_cover)
  dl_4_cover <- dl_4_cover[-grep('hold', names(dl_4_cover))]
  names(dl_4_cover) <- paste0( 'cover_', names(dl_4_cover))
  
  if(vr == 'growth'){ 
    dat <- dat[complete.cases(dat), ]
    dat <- split_df(dat, ... )
    dl <- make_dl(dat)
  }else if(vr == 'survival'){ 
    dat <- split_df(dat, ... )
    dl  <- make_dl(dat)
  }
  
  return( c(dl, dl_4_cover))
}


process_recruitment_data <- function(dat, formX, formC, formZ, center = T, ... ){ 
  
  C <- model.matrix(formC, dat)
  dat$C <- scale(C)
  dat$Group <- factor(dat$gid)
  
  dat$X <- model.matrix(formX, data = dat)
  dat$Z <- model.matrix(formZ, data = dat)
  
  dat$g <- factor(dat$yid)
  
  dat <- split_df(dat, hold )
  
  dl <- make_dl(dat)
  
  return(dl)
}


check_for_compiled_model <- function(vr, model_file){ 
  compiled_model <- dir('data/temp_data/', paste0(vr, '_model_compiled.RDS'), full.names = T)  
  while(length(compiled_model) == 0){
    rt <- stanc(model_file)
    stan_mod <- stan_model(stanc_ret = rt, verbose = FALSE, auto_write = F)
    saveRDS(stan_mod, file.path('data/temp_data', paste0(vr, '_model_compiled.RDS')))
    compiled_model <- dir('data/temp_data', paste0(vr, '_model_compiled.RDS'), full.names = T)
  }
  return( readRDS(compiled_model) )
}

drop_init_years <- function(my_inits, hold){ 
  
  lapply( my_inits, function(x){x$u_raw = x$u_raw[, -hold]; return(x) } )
}

drop_init_covariates <- function(my_inits, K){
  lapply( my_inits, function(x) { x$beta <- x$beta[1:K]; return(x) })  
}

find_dv_trans <- function(x){ 
  ss <-  get_sampler_params(x) 
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  return(dv)
}

