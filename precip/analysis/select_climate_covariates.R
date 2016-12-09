# 
library(stringr)
rm(list = ls())

cor_files <- dir('output', 'correlations.csv', full.names = T)

out <- list( NA )

for(i in 1:length(cor_files)){

  temp <- read.csv(cor_files[i])
  
  parameter <- as.character( temp$X)
  size_par <- which( temp$size_pval < 0.1)
  
  if(length(size_par) != 0 ){
    ifx <- paste( parameter[ size_par ], 'logarea.t0', sep = 'x')
    use_vars <- c( parameter, ifx)
  }else{ 
    use_vars <- parameter
  }
  
  spp <- str_extract( pattern =  '[A-Z]{4}', cor_files[i])
  vr  <- str_extract( pattern = '(growth)|(recruitment)|(survival)', cor_files[i])
  
  out[[i]] <- data.frame( species = spp, vital_rate = vr, covars = paste( use_vars, collapse = ','))
}

out <- do.call(rbind, out)

write.csv(out, 'output/selected_climate_covariates.csv')

# check for highly correlated covariates 
clim <- readRDS('data/temp_data/all_clim_covs.RDS')
clim <- clim[complete.cases(clim), ]
correlations <- list()
i = 3
for ( i in 1:nrow( out )  ) { 
  covars     <- strsplit( as.character( out$covars[i] ) , ',')[[1]]
  covars <- covars[ -grep('x', covars)]
  correlations[[i]] <- cor(clim[ , covars, drop = F] )
}

corelations 



