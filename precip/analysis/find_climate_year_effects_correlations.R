rm(list = ls())
library(ggplot2)
library(stringr)

ye_files <- dir('output', 'year_effects_table.csv', full.names = T)
dat_files <- dir('data/temp_data', 'data_lists_for_stan.RDS', full.names = T, recursive = T)

Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')
Cdat <- subset(Cdat, Period == 'Historical' & Treatment == 'Control')

for( i in 1:length(ye_files)){ 
  
  spp <- str_extract (ye_files[i], '[A-Z]{4}')
  vr  <- str_extract (ye_files[i], c('growth', 'recruitment', 'survival'))
  vr  <- vr[ !is.na(vr)]
  
  ye <- read.csv(ye_files[i])
  
  ye$parameter  <- str_match(ye$X, '(^[a-z]+1?)\\[')[, 2]
  
  ye <-  subset( ye , parameter %in% c('a', 'b1'))  
  
  dl <- readRDS( dat_files [ str_detect(dat_files, vr) ] )[[spp]]

  years <- unique(dl$year)

  ye_df <- data.frame( year = years, parameter = ye$parameter, year_effect = ye$mean )
  
  ye_df <- ye_df[order(ye_df$parameter, ye_df$year), ]

  pdf( paste0( 'figures/year_effects_for', spp, '_', vr, '.pdf' ))
  print( 
    ggplot( ye_df, aes( x = year, y = year_effect, color = parameter , linetype = parameter)) + 
      geom_point() + 
      geom_line() + 
      ggtitle(paste0('Year effects for ', spp, ' ', vr))
  )
  
  ggplot( ye_df, aes( x = year, y = year))
  
  dev.off()
  
  ye_df <- ye_df %>% spread(parameter, year_effect)
  
  ye_df <- merge(ye_df, Cdat, by = 'year')
  
  clim <- ye_df[ ,grep('^T\\.|^VWC\\.', names(ye_df)) ]
  
  if ( vr == 'recruitment'){ 
    
    out <- data.frame(matrix(NA, nrow = ncol(clim), ncol = 2))
    names(out) <- c('Intercept', 'Int_pval')
    row.names(out) <- names(clim)
    
    for(j in 1:ncol(clim)){ 
      
      temp <- cor.test(clim[, j], ye_df$a)
      out [ j, 1] <- temp$estimate
      out [ j, 2] <- temp$p.value
  
    }
    
    out <- subset(out, Int_pval < 0.1)
    out <-  out[order(out$Intercept), ]
    write.csv(out, paste0( 'output/', spp, '_', vr, '_correlations.csv'))      
  }else{ 
  
    out <- data.frame(matrix(NA, nrow = ncol(clim), ncol = 4))
    
    names(out) <- c('Intercept', 'Int_pval', 'size', 'size_pval')
    row.names(out) <- names(clim)
    
    for(j in 1:ncol(clim)){ 
      
      temp <- cor.test(clim[, j], ye_df$a)
      out [ j, 1] <- temp$estimate
      out [ j, 2] <- temp$p.value
      
      temp <- cor.test(clim[, j], ye_df$b)
      out [ j, 3] <- temp$estimate
      out [ j, 4] <- temp$p.value
    }
    
    out <- subset(out, Int_pval < 0.1 | size_pval < 0.1)
    out <-  out[order(out$Intercept), ]
    write.csv(out, paste0( 'output/', spp, '_', vr, '_correlations.csv'))  
  } 
} 
  