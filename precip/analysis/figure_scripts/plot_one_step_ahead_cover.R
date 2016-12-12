rm(list =ls())
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)

# functions ---------------------------------------------------------------------------- 

simulate_recruit_size <- function( recruits , spp ) { 
  
  dat <- read.csv(paste0('~/driversdata/data/idaho/speciesData/', as.character(spp), '/recSize.csv'))
  
  out <- recruits 
  
  for(i in 1:length(recruits)){ 
    # draw size distribution of n recruits in each plot i 
    # sum up total size to get total cover of recruits in plot i  
    n <- round(recruits[i])
    out[i] <- sum ( sample(dat$area, n, replace = T) )  
  }
  return(  out  )
}
  
#
load('analysis/figure_scripts/my_plotting_theme.Rdata') 

years <- expand.grid(year = 1925:2017, Treatment = c('Control', 'Drought', 'Irrigation'), stat = c('observed', 'predicted'))
years$Period[ years$year > 2006 ] <- 'Modern'
years$Period[ years$year <= 1960 ] <- 'Historical'

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

ylims <- list( c(0,35), c(0,7.5), c(0,7.5), c(0,7.5))
i = 1


for( i in 1:4) {  
  spp   <- species_list[i]  
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[i]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[i]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[i]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_climate_fit.RDS'), full.names = TRUE)
  
  s <- rstan::extract(readRDS(m[3]), c('muhat2'))$muhat2
  g <- rstan::extract(readRDS(m[1]), c('muhat3'))$muhat3
  r <- rstan::extract(readRDS(m[2]), c('lambda_pred2'))$lambda_pred2
  
  df <- read.csv(paste0( 'data/temp_data/', spp, '_growth_and_survival_cleaned_dataframe.csv'))
  rdf <- read.csv(paste0( 'data/temp_data/', spp, '_recruitment_cleaned_dataframe.csv'))
  
  rec_area <- simulate_recruit_size(r, spp )
  g <- exp( g ) 
  g[ g > 10000 ]  <- 10000 # assign plants with greater than 10000 cover to 10000 (can't be bigger than plot)
  k <- s*g 
  
  df <- data.frame(df, t(k))
  rdf <- data.frame(rdf, t(rec_area))
  
  cover1 <-   
    df %>% 
    gather( simulation, size , matches('^X[0-9]')) %>% 
    group_by( year, Treatment, simulation, quad ) %>% 
    summarize( total_size = sum( size )) 

  cover2 <- 
    rdf %>% 
    gather( simulation, rec_area, matches('^X[0-9]')) %>% 
    group_by( year, Treatment, simulation) %>% 
    summarise( rec_area = mean(rec_area))
  
  cover <- 
    merge( cover2, cover1, all.y = TRUE) %>%  # use all recruitment years  
    mutate( total_size = ifelse(is.na(total_size), 0, total_size)) %>%  # years without observed plants get 0
    mutate( area = total_size + rec_area ) %>% 
    group_by( year , Treatment, simulation ) %>% 
    summarise(cover = 100*( mean(area)/10000) ) # percent cover is total area / 10000 cm * 100

  pred_cover <- 
    cover %>% 
    ungroup %>% 
    mutate( year = year + 1 ) %>% # cover predictions the "Y's" are for the next year 
    group_by(Treatment, year ) %>% 
    summarise( predicted = mean(cover), 
               ucl = quantile( cover, 0.75), 
               lcl = quantile(cover , 0.25))
  
  # get last year of cover which is not in the survival dataframe ------------------------------------- #   
  oldCover <- read.csv(paste0( '~/driversdata/data/idaho/speciesData/', spp, '/quadratCover.csv'))
  newCover <- read.csv(paste0( '~/driversdata/data/idaho_modern/speciesData/', spp, '/quadratCover.csv'))
  oldCover$Period <- "Historical"
  newCover$Period <- "Modern"
  last_cover  <- rbind(oldCover, newCover)
  quad <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')
  
  last_cover <-
    last_cover %>%
    left_join(quad) %>%
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>%
    mutate( year = ifelse(year < 100, year + 1900, year)) %>%
    group_by(Period, Treatment, year) %>%
    summarise( true_cov = 100*(mean(totCover )/10000)) 
  
  # ------------------------------------------------------------------------------------------------------ #
  df <- read.csv('data/temp_data/ARTR_growth_and_survival_cleaned_dataframe.csv')
  
  df$true_cov <- 100*(exp(df$logarea.t0)/10000)
  
  trueCov <- 
    df %>% 
    group_by( Period, year, Treatment, quad ) %>% 
    summarise( totCover = sum( true_cov ) ) %>% 
    group_by( Period, year, Treatment) %>% 
    summarise( observed  = mean(totCover ))
  
  obs_cover <- merge( last_cover, trueCov , by = c('Period', 'year', 'Treatment'), all.x = T)
  
  obs_cover$observed <- ifelse( is.na(obs_cover$observed), obs_cover$true_cov, obs_cover$observed) # fill in data from end years in dataframe
  
  rm(sdl, rdl)
  
  plot_cover <-
    merge( obs_cover, pred_cover, all.x = TRUE) %>%
    mutate( predicted = ifelse(Treatment != 'Control' & year == 2011, observed, predicted )) %>%
    mutate( predicted = ifelse(year == 2007, observed, predicted )) %>%
    mutate( predicted = ifelse(year == 1929, observed, predicted )) %>%
    mutate( predicted = ifelse(year == 1945, observed, predicted )) %>%
    mutate( predicted = ifelse(year == 1949, observed, predicted )) %>%
    mutate( predicted = ifelse(year == 1954, observed, predicted )) %>%
    mutate( ucl = ifelse(Treatment != 'Control' & year == 2011, observed, ucl )) %>%
    mutate( ucl = ifelse(year == 2007, observed, ucl )) %>%
    mutate( ucl = ifelse(year == 1929, observed, ucl )) %>%
    mutate( ucl = ifelse(year == 1945, observed, ucl )) %>%
    mutate( ucl = ifelse(year == 1949, observed, ucl )) %>%
    mutate( ucl = ifelse(year == 1954, observed, ucl )) %>%
    mutate( lcl = ifelse(Treatment != 'Control' & year == 2011, observed, lcl )) %>%
    mutate( lcl = ifelse(year == 2007, observed, lcl )) %>%
    mutate( lcl = ifelse(year == 1929, observed, lcl )) %>%
    mutate( lcl = ifelse(year == 1945, observed, lcl )) %>%
    mutate( lcl = ifelse(year == 1949, observed, lcl )) %>%
    mutate( lcl = ifelse(year == 1954, observed, lcl )) %>%
    gather( stat, val,  predicted, observed) 
  
  plot_cover <- merge( years, plot_cover, all.x = TRUE)
  
  pdf( paste0( 'figures/predictions/', spp  , '_predicted_cover.pdf' ), height = 8, width = 8) 
  
  print( 
    ggplot( subset( plot_cover, Period == "Historical"), aes( x = year, y =  val, fill = Treatment, color = Treatment, linetype = stat, ymax = ucl, ymin = lcl)) +
      geom_line() +
      #geom_ribbon(aes(ymax = ucl, ymin = lcl, color = NA), alpha = 0.1) + 
      scale_color_manual(values = my_colors ) + 
      ylim( ylims[[i]]) + 
      ylab( 'cover (%)' ) +  
      scale_x_continuous(breaks = seq(1925,1960,5)) + 
      my_theme + 
      guides( linetype=guide_legend(title=NULL) ) + 
      ggtitle(paste('Historical cover of', spp )) 
  )
  
  print( 
    ggplot( subset( plot_cover, Period == "Modern"), aes( x = year, y =  val, fill = Treatment, color = Treatment, linetype = stat, ymax = ucl, ymin = lcl)) +
      geom_line() +
      #geom_ribbon(aes(ymax = ucl, ymin = lcl, color = NA), alpha = 0.1) + 
      scale_color_manual(values = my_colors ) + 
      ylim( ylims[[i]] ) + 
      ylab( 'cover (%)') + 
      scale_x_continuous(breaks = c(2007:2016)) + 
      my_theme + 
      guides( linetype=guide_legend(title=NULL) ) + 
      ggtitle(paste('Modern cover of', spp ))
  )
  
  dev.off()
} 

saveRDS(size , 'output/predicted_size.RDS')
saveRDS(survives, 'output/predicted_survival.RDS')
