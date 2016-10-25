tweak_inits <- function(inits, devs ){ 
  # adds random deviation to list of inital values 
  inits + runif(length(inits), 0, 0.5)
}

run_stan_model <- function(do_spp, do_vr, do_model, do_lambda, do_prior_sd, pars = NULL, nchains = 1, niter = 1, predict = FALSE) { 
  
  require(rstan)
  
  data_path <- file.path(getwd(), 'data/temp_data')
  model_path <- file.path(getwd(), 'analysis')
  
  if( predict )  { 
    output_path <- file.path(getwd(),  'output/stan_fits/predictions')
  }else{ 
    output_path <- file.path( getwd(), 'output/stan_fits')
  }
  
  # read in data and models ---------------------------------------------------------------------------------------#
  do_stan_file <- do_model

  data_file <- dir(data_path, pattern = paste0(do_vr, '_data_lists_for_stan.RDS'), full.names = TRUE )
  
  if ( predict ){ 
    initial_fit <- dir(output_path, pattern = paste(do_spp, do_vr, do_model, '[0-9]+', 0, 'predict.RDS', sep = '_'), full.names = TRUE) # check pre-fit models
  }else{ 
    initial_fit <- dir(output_path, pattern = paste(do_spp, do_vr, do_model, '[0-9]+', 0, '.RDS', sep = '_'), full.names = TRUE) # check pre-fit models
  }
  
  data_list <- readRDS( data_file )[[do_spp]]
  
  if(length(data_list) == 0) { stop('No data!!!')}

  if(predict) { 
    models <- dir(file.path(model_path, do_vr), '[0-9]_predict.stan', full.names = TRUE)
  } else { 
    models <- dir(file.path(model_path, do_vr), '[0-9].stan', full.names = TRUE)
  }
  
  
  # -- select model and initial vals -----------------------------------------------------------------------------# 
  
  m <- models [do_stan_file]
  
  # modify data list for model ------------------------------------------------------------------------------------#
  
  data_list$tau_beta <- do_prior_sd  # modify prior for regularization 
  
  # #
  # # select only intraspecific w's for models with intraspecific crowding only 
  # #
  
  if( !do_vr == 'recruitment') {
      if( do_model < 3 ) {  ## intraspecific effects only 
        
        data_list$W <- data_list$W[ ,  data_list$spp ]
        data_list$Whold <- data_list$Whold[ , data_list$spp ]
        data_list$W2 <- data_list$W2[ , data_list$spp ]
        data_list$W3 <- data_list$W3[ , data_list$spp ]
      }
  }
  
  # get initial values for proper species and model --------------------------------------------------------------# 
  
  init_file <- dir(data_path, pattern = paste0(do_vr, '_init_vals.RDS'), full.names = TRUE)
  init_vals <- readRDS(init_file)[[do_spp]]
  temp_inits <- init_vals[[do_model]]
  temp_inits <- rep(list(temp_inits), max(1, nchains) )
  devs <- seq(-0.5,0.5,length.out = nchains)
  temp_inits <- lapply(temp_inits, lapply, tweak_inits)
  
  # -- run stan ---------------------------------------------------------------------------------------------------#  
  
  save_file <- file.path( output_path, paste(do_spp, do_vr, do_model, do_lambda, nchains, sep = '_'))
  
  print(paste('data list has', length(data_list), 'variables'))
  
  if ( length(initial_fit) > 0 ) { 
    print(paste('initial fit being used', initial_fit ))
    initial_fit <- readRDS(head( initial_fit ) )
    temp_fit <- stan (fit = initial_fit, model_name = basename(save_file), data = data_list, chains = nchains, iter = niter, pars = pars, init = temp_inits, cores = max(1, nchains))
  }else{   
    temp_fit <- stan (file = m, model_name = basename(save_file), data = data_list, chains = nchains, iter = niter , pars = pars, init = temp_inits, cores = max(1, nchains))
  }
  
  # -- output -----------------------------------------------------------------------------------------------------#
  
  return(temp_fit)
  
}