# 
# 
# get WAIC for missing runs  
# 
# 

rm(list = ls() ) 

library(rstan)

args <- commandArgs(trailingOnly=TRUE)
#args <- c('~/Documents/ExperimentTests/precip/' , '1', '1')

# test if there is at least one argument: if not, return an error
if (length(args) != 3){ 
  stop('####### Incorrect number of arguments supplied ####### \n
       ####### Arguments required:
       #######  working directory 
       #######  line number : 1 - total combination of models in
       #######  chains = 4 hardcoded  
       #######  niter = 2000 hardcoded')
  
}else if (length(args) == 3){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  models <- read.csv(paste0('output/missing_runs_oos', args[2], '.csv'))

  do_line <- as.numeric(eval(parse(text = args[3])))
  
  # nchains <- as.numeric(eval(parse (text = strsplit( args[3], ' '))))
  # niter <- as.numeric(eval(parse (text = strsplit( args[4], ' '))))
  
}

nchains <- 4 
niter <- 2000


source( 'analysis/run_stan_fxns.R')
source( 'analysis/waic_fxns.R')

if ( do_line <= nrow(models)) { 
  
  line <- models[do_line, ]
  
  species <- line$species 
  vital_rate <- line$vital_rate
  model <- line$model
  prior <- line$sd
  lambda <- line$lambda
  nlambda <- line$nlambda
  year_oos <- as.character( line$year_oos)
  
}else{ stop('line number is greater than number of models')}

output_path <- file.path(getwd(),  'output/stan_fits/WAIC_scores/oos_scores')

save_file <- file.path( output_path, paste(species, vital_rate, model, lambda, nchains, year_oos, 'WAIC.csv', sep = '_'))

set.seed(9)
temp_fit <- run_stan_model_oos(species, vital_rate, model, lambda , prior, year_oos = year_oos, nchains = nchains, niter = niter, pars = c('log_lik', 'log_lik2'), predict = FALSE)

waic_df <- rbind( waic(temp_fit, llname = 'log_lik'), waic(temp_fit, llname = 'log_lik2'))
waic_df$type <- c('in_sample', 'out_of_sample')

waic_df$fn <- basename(save_file)

write.csv(waic_df, file = save_file, row.names = FALSE)

