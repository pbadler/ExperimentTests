
# 
library(lme4)
library(stringr)

mydf <- dir( 'data/temp_data/', '[A-Z]{4}_growth_cleaned_dataframe.RDS' , full.names = T)

growth_f <- as.formula(Y ~ X + (X|year) + Group + W.ARTR + W.HECO + W.POSE + W.PSSP + 
                         T.sp.0 + T.sp.1 + VWC.sp.0 + VWC.sp.1 + VWC.f.1 + VWC.f.0 + T.f.1 + T.f.0 )


fit <- list(NA)
x   <- list(NA)
cover_df <- list(NA)
for(i in 1:length(mydf) ){  
  x[[i]] <- readRDS(mydf[i])
  species <- str_extract(basename(mydf[i]), '[A-Z]{4}')
  
  fit[[i]] <- lmer( data = subset(x[[i]], Period == 'Historical'), growth_f)
  
  x[[i]]$predicted <- predict( fit[[i]], newdata = x[[i]], re.form = NA)
  x[[i]]$diff <- x[[i]]$predicted - x[[i]]$Y
  
  pdf( paste0( 'figures/lmer_predictions_', species, '_', 'growth', '.pdf' ))
    par(mfrow = c(2,1))
    temp <- subset(x[[i]], Period == 'Modern')
    plot(temp$predicted, temp$Y, xlab = 'Predicted', ylab = 'Observed', main = paste0(species, ' growth')) 
    abline( 0, 1, lwd = 2, col = 'red')
  
    hist(temp$diff, xlab = 'Predicted - Observed', main = '') 
    abline(v = 0 , col = 'red', lwd = 5)
  dev.off()
  
  
  MSE_by_year <- x[[i]] %>% 
    group_by(Period, Treatment, year) %>% 
    summarise( mse = mean(diff^2), msd = mean(diff))
  
  MSE_by_treatment <- x[[i]] %>% 
    group_by(Period, Treatment) %>% 
    summarise( mse = mean(diff^2), msd = mean(diff))
  
  MSE_by_period <- x[[i]] %>% 
    group_by(Period) %>% 
    summarise( mse = mean(diff^2), msd = mean(diff))
  
  # cover predictions 
  
  cover_df[[i]] <- readRDS(paste0( 'data/temp_data/', species, '_survival_cleaned_dataframe.RDS' ))
  
  cover_df[[i]]$predicted <-  predict( fit[[i]], cover_df[[i]], re.form = NA)
  
}

saveRDS(x, 'output/lmer_predictions/lmer_growth.RDS')
saveRDS(cover_df, 'output/lmer_predictions/lmer_growth_cover.RDS')


lapply( fit, summary)
