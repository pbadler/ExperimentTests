rm(list  = ls() )

library(rstan)
library(tidyverse)

# Extract growth and survival year effects first 

sp_vr <- expand.grid( spp = c('ARTR', 'HECO', 'POSE', 'PSSP'), vr = c('growth', 'survival'))

out <- list() 

for(i in 1:nrow(sp_vr)){ 
    
  sp <- sp_vr$spp[i]
  vr <- sp_vr$vr[i]
  
  fit <- readRDS(paste0( '~/Desktop/', sp, '_', vr, '_fit1.RDS'))
  dat <- readRDS(paste0( '~/Desktop/', sp, '_', vr, '_data.RDS'))

  year_effects <- summary(fit, 'u')$summary

  nr <- nrow(year_effects)

  if(vr != 'recruitment' ) { 
    int <- year_effects[seq(1,nr,2), ]
    slope <- year_effects[1 + seq(1,nr,2), ]
  
    int_df <- data.frame(type = 'intercept', year = dat$years, int[, c(1,4,8)])
    slope_df <- data.frame(type = 'slope', year = dat$years, int[, c(1,4,8)])

    year_effects <- rbind(int_df, slope_df)
    
  }else if (vr == 'recruitment' ) { 
  
    year_effects <- data.frame(type = 'intercept', year = dat$years, year_effects[, c(1,4,8)])
  }
  
  year_effects$spp <- sp 
  year_effects$vr  <- vr 
  
  out[[i]] <- year_effects
}

write_csv(path = 'output/year_effects.csv' , do.call( rbind, out ) )
