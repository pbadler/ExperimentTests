# 
# 
# get WAIC simple 
# 
# 

rm(list = ls() ) 

library(rstan)

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 4){ 
  stop('####### Incorrect number of arguments supplied ####### \n
       ####### Arguments required:
       #######  working directory 
       #######  line number : 1 - total combination of models in 
       #######  chains: 1 - 4 
       #######  niter:  number of iterations to run per chain')
  
}else if (length(args) == 4){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  do_line <- as.numeric(eval(parse(text = args[2])))
  
  nchains <- as.numeric(eval(parse (text = strsplit( args[3], ' '))))
  niter <- as.numeric(eval(parse (text = strsplit( args[4], ' '))))
  
}

models <- read.csv('data/temp_data/model_table.csv')
source( 'analysis/run_stan_fxns.R')
source( 'analysis/waic_fxns.R')


if ( do_line < nrow(models)) { 
  
  line <- models[do_line, ]
  
  species <- line$species 
  vital_rate <- line$vital_rate
  model <- line$model
  prior <- line$prior
  nlambda <- line$nlambda

}else{ stop('line number is greater than number of models')}

output_path <- file.path(getwd(),  'output/stan_fits/WAIC_scores/')

save_file <- file.path( output_path, paste(species, vital_rate, model, prior, nchains, 'WAIC.csv', sep = '_'))
temp_fit <- run_stan_model(species, vital_rate, model, prior, nchains = nchains, niter = niter, pars = 'log_lik', predict = FALSE, nlambda = nlambda )

write.csv(waic(temp_fit), file = save_file, row.names = FALSE)

