# 
library(stringr)
rm(list = ls())

cor_files <- dir('output', 'correlations.csv', full.names = T)
cor_files <- cor_files [ -grep('all', cor_files)] # keep only the selected ones 

out <- list( NA )


for(i in 1:length(cor_files)){
  
  spp <- str_extract( pattern =  '[A-Z]{4}', cor_files[i])
  vr  <- str_extract( pattern = '(growth)|(recruitment)|(survival)', cor_files[i])
  
  temp <- read.csv(cor_files[i])
  
  VWC_parameter <- as.character( temp$var [ temp$vartype == 'VWC'] )
  T_parameter <- as.character( temp$var[ temp$vartype == 'T'])
  
  if ( vr == 'growth' & spp %in% c('ARTR', 'POSE')) {  # use the size effect interaction 
    ifx <- paste( VWC_parameter, 'logarea.t0', sep = 'x')
    
    VWC_parameter <- c(VWC_parameter, ifx )
  }

  if( vr == 'survival' & spp == 'POSE' ) { # use 
    ifx <- paste(VWC_parameter, 'logarea.t0', sep = 'x')
    VWC_parameter <- c(VWC_parameter, ifx)
  }  
  
  use_vars <- c(T_parameter, VWC_parameter)
  
  out[[i]] <- data.frame( species = spp, vital_rate = vr, covars = paste( use_vars, collapse = ','))
}

out <- do.call(rbind, out)

write.csv(out, 'output/selected_climate_covariates.csv')

# check for highly correlated covariates 
# clim <- readRDS('data/temp_data/all_clim_covs.RDS')
# clim <- clim[complete.cases(clim), ]
# correlations <- list()
# i = 3
# for ( i in 1:nrow( out )  ) { 
#   covars     <- strsplit( as.character( out$covars[i] ) , ',')[[1]]
#   covars <- covars[ -grep('x', covars)]
#   correlations[[i]] <- cor(clim[ , covars, drop = F] )
# }
# 
# all_cors <- lapply(cor_files,  read.csv)
# 
# names(all_cors) <- basename( cor_files )
# 
# recruitment_cors <- do.call( rbind, all_cors[ grep( 'recruitment', names(all_cors)) ] )
# 
# labels <- do.call( rbind, str_split(row.names(recruitment_cors), '_'))[, 1:2]
# recruitment_cors$species <- labels[, 1]
# recruitment_cors$vital_rate <- labels[, 2]
# 
# gs_cors <- do.call( rbind, all_cors[ -grep( 'recruitment', names(all_cors)) ] )
# 
# labels <- do.call( rbind, str_split(row.names(gs_cors), '_'))[, 1:2]
# gs_cors$species <- labels[, 1]
# gs_cors$vital_rate <- labels[, 2]
# 
# recruitment_cors
# gs_cors
# 
# write.csv(file = 'output/year_effects_correlations_recruitment.csv', recruitment_cors, row.names = F)
# write.csv(file = 'output/year_effects_correlations_growth_survival.csv', gs_cors, row.names = F)
# 


