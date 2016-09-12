###########################################################################
#
##  Functions to automate running stan models 
##  
##  Runs separate models for each species and vital rate. Only works
##  for models that do not involve climate effects, 
##  i.e. do not require regularization. 
##
##  Saved output:
##    Save full model output from each model  
##
# 
##########################################################################


tweak_inits <- function(inits){ 
  
  inits + runif(length(inits), -0.001, 0.001)
  
}

run_stan_model <- function(do_spp, do_vr, do_model, nchains, niter ) { 
  
  require(rstan)
  
  data_path <- file.path(getwd(), 'data/temp_data')
  model_path <- file.path(getwd(), 'analysis')
  output_path <- file.path(getwd(),  'output/stan_fits')
  
  # read in data and models ---------------------------------------------------------------------------------------#
  
  data_file <- dir(data_path, pattern = paste0(do_vr, '_data_lists_for_stan.RDS'), full.names = TRUE )
  init_file <- dir(data_path, pattern = paste0(do_vr, '_init_vals.RDS'), full.names = TRUE)

  data_list <- readRDS( data_file )[do_spp][[1]]
  init_vals <- readRDS(init_file)[do_spp][[1]]
  
  models <- dir(file.path(model_path, do_vr), '[0-9].stan', full.names = TRUE)
  
  # -- select model and initial vals -----------------------------------------------------------------------------# 
  

  m <- models [do_model]
  
  temp_inits <- init_vals[[do_model]]
  
  # modify data list for model ------------------------------------------------------------------------------------#
  
  #
  # select only intraspecific w's for models with intraspecific crowding only 
  #
  
  if( 'w' %in% names( temp_inits )  ) { 
    if( length(temp_inits$w) == 1) { 
      w_names <- colnames(data_list$W)
      data_list$W <- data_list$W[ ,  grep(pattern = do_spp, w_names) ] 
      data_list$Whold <- data_list$Whold[ , grep(pattern = do_spp, w_names) ]
      data_list$W_covs <- 1 
    } 
  }
  
  # -- inital values ----------------------------------------------------------------------------------------------#
  
  temp_inits <- rep(list(temp_inits), nchains )
  temp_inits <- lapply(temp_inits, lapply, tweak_inits)
  
  # -- run stan ---------------------------------------------------------------------------------------------------#  
  
  save_file <- file.path( output_path, paste(do_spp, do_vr, do_model, sep = '_'))
  
  temp_fit <- stan (file = m, model_name = basename(save_file), data = data_list, chains = nchains, init = temp_inits, iter = niter , cores = max(1, nchains))
  
  # -- output -----------------------------------------------------------------------------------------------------#
  
  saveRDS(temp_fit, file = paste0( save_file, '.RDS'))
  
}
