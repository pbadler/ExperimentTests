library(lme4)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

Cdat <- read.csv('~/driversdata/data/idaho/climateData/olderAdlerFiles/Climate.csv')
VWCdat <- readRDS('data/temp_data/all_clim_covs.RDS')
VWCdat <- subset( VWCdat, Period == 'Historical' & Treatment == 'Control')

pairs(VWCdat[, grep('^T\\.|VWC\\.', names(VWCdat))])

Cdf <- merge( Cdat, VWCdat)

pairs(Cdf[, c('pptLag', 'ppt1', 'TmeanSpr1', 'ppt2', 'TmeanSpr2', 'VWC.sp.l', 'VWC.sp.0', 'VWC.sp.1', 'VWC.su.0')])
pairs(Cdf[, c('ppt1', 'VWC.sp.0', 'P.f.w.sp.0')])

Gfiles <- dir('data/temp_data', '*_growth.RDS', full.names = T)

for(i in 1:length(Gfiles)){
  
  G  <-  readRDS(Gfiles[i])
  #df <- merge(G, Cdat, by = 'year')
  df <- merge(G, VWCdat, by = 'year')

  
  myr <- lmer(data = df, 
              logarea.t1 ~ 
                logarea.t0 + 
                W.ARTR + 
                W.HECO + 
                W.POSE + 
                W.PSSP + 
                (logarea.t0|year) + 
                (1|Group))
  
  #year_effects <- cbind( unique(df[, grep('^ppt', names(df))]), ranef(myr)$year)
  VWC_year_effects <- cbind( unique(df[, grep('^VWC', names(df))]), ranef(myr)$year)
  T_year_effects <- cbind( unique(df[, grep('^T\\.', names(df))]), ranef(myr)$year)
  
  spp <-  str_extract(Gfiles[i], '[A-Z]{4}')  
  
  pdf(paste0('figures/growth_correlation_year_effects' , spp , '.pdf') , width = 8, height = 4)
  par(mfrow=c(1,2))

  for(j in 1:ncol(VWC_year_effects)){
    
    plot(VWC_year_effects[,j ], VWC_year_effects$`(Intercept)`, xlab = names(VWC_year_effects)[j])
    plot(VWC_year_effects[,j ], VWC_year_effects$logarea.t0, xlab = names(VWC_year_effects)[j])
    
  }
  
  for(j in 1:ncol(T_year_effects)){
    
    plot(T_year_effects[,j ], T_year_effects$`(Intercept)`, xlab = names(T_year_effects)[j])
    plot(T_year_effects[,j ], T_year_effects$logarea.t0, xlab  = names(T_year_effects)[j])
    
  }
  
  
  dev.off()

  ye <- merge (T_year_effects, VWC_year_effects)
  
  out <- data.frame(matrix(NA, nrow = ncol(ye), ncol = 4))
  
  names(out) <- c('(Intercept)', 'Int_pval', 'logarea', 'logarea_pval')
  row.names(out) <- names(ye)

  for(j in 1:ncol(ye)){ 
    
    Rest <- cor.test(ye[, j], ye$`(Intercept)`)
    Rest$estimate
    Rest$p.value
    
    out [ j, 1] <- Rest$estimate
    out [ j, 2] <- Rest$p.value
    
    Rest <- cor.test(ye[, j], ye$logarea.t0)
    out [ j, 3] <- Rest$estimate
    out [ j, 4] <- Rest$p.value
  }
  
  out <- subset(out, Int_pval < 0.1 | logarea_pval < 0.1)
  out <-  out[order(out$`(Intercept)`), ]
  
  write.csv(out, paste0( 'output/', spp, '_growth_correlations.csv'))
}

