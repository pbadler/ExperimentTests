rm(list =ls())
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)
library(rstan)

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
  return(out )
}


#
setwd('~/Documents/ExperimentTests/precip/')

my_colors <- c('#1b9e77', '#d95f02', '#7570b3')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

ylims <- list( c(0,40), c(0,7.5), c(0,7.5), c(0,7.5))
i = 2

for( i in 1:4) {  
  spp   <- species_list[i]  
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[i]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[i]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[i]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_climate_fit.RDS'), full.names = TRUE)
  
  temp_fit <- readRDS(m[3])
  s1 <- summary(temp_fit, 'mu')$summary[, 1]
  s2 <- summary(temp_fit, 'muhat')$summary[, 1]
  s <- c(s1, s2)
  
  g <- summary(readRDS(m[1]), 'muhat3')$summary[,1]
  g[ exp(g)/10000 > 1 ] <-  log(10000)  # assign plants larger than the entire plot to the plot size 
  
  r1 <- summary(readRDS(m[2]), 'lambda')$summary[,1]
  r2 <- summary(readRDS(m[2]), 'lambda_pred')$summary[,1]
  r <- c( r1, r2 )

  rsize <- simulate_recruit_size( r, spp = spp)
  rsize[ (rsize > 10000 ) ] <- 10000 # assign cover over 10000 (100%) to 10000
  
  nsims <- dim(s)[1]
  nobs  <- dim(s)[2]
  
  # combine survival, growth and recruitment ---------------------------------------------------------- # 
  g <- 100*( exp( g ) /10000 ) # convert to % cover 
  size <- s*g                  # survival by size 
  rsize <- 100*(rsize)/10000   # convert recruits to % cover 
  
  size <- data.frame( Period = gdl$Period3, year = gdl$year3, Treatment = gdl$treat3 , quad = gdl$quad3, Group = gdl$gid3 , size = size)
  rsize  <- data.frame(Period = rdl$Period2, year = rdl$year2, Treatment = rdl$treat2 , quad = rdl$quad2, Group = rdl$gid2, rsize = rsize)
  
  cover1 <-  
    size %>% 
    group_by( Period, Group, Treatment, year, quad) %>% 
    summarise( size_cover = sum( size ))
  
  cover <- merge(rsize, cover1, all.x = T)
  cover$size_cover [ is.na( cover$size_cover) ] <- 0 # assign size cover to zero when no plants were present 
  cover$pred_cover <- cover$size_cover + cover$rsize
  cover$pred_cover[ cover$pred_cover > 100 ] # total cover can be greater than 100 because of overlapping canopies
  cover$Treatment <- factor(cover$Treatment, labels = c('Control', 'Drought', 'Irrigation'))
  
  # get last year of cover which is not in the survival dataframe ------------------------------------- #   
  oldCover <- read.csv(paste0( '~/driversdata/data/idaho/speciesData/', spp, '/quadratCover.csv'))
  newCover <- read.csv(paste0( '~/driversdata/data/idaho_modern/speciesData/', spp, '/quadratCover.csv'))
  oldCover$Period <- 'Historical'
  newCover$Period <- 'Modern'
  obs_cover  <- rbind(oldCover, newCover)
  quad <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')
  
  obs_cover$obs_cover <- (100*obs_cover$totCover)/10000 # convert to percent cover 
  obs_cover <- merge( obs_cover, quad)
    
  obs_cover <-
    obs_cover %>%
    left_join(quad) %>%
    filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>%
    mutate( year = ifelse(year < 100, year + 1900, year)) %>% 
    mutate( quad = str_extract(quad, '[0-9]+'))
  
  test <-  merge(obs_cover, cover [, c('year', 'pred_cover', 'Treatment', 'quad')], by = c('year' , 'quad', 'Treatment') )

  ggplot( test, aes( x = pred_cover, y = obs_cover) ) + geom_point() + facet_grid( . ~ Period ) + geom_abline(aes( intercept = 0 , slope = 1))
  
  ggplot( test, aes( x = pred_cover, y = obs_cover) ) + geom_point() + 
    facet_grid( Period ~ Group ) + 
    geom_abline(aes( intercept = 0 , slope = 1))
  
  avg_cover <- 
    test %>% 
    gather( type , cover, pred_cover, obs_cover ) %>% 
    group_by( Period, year , Group, Treatment, type ) %>% 
    summarise( cover = mean(cover)) %>% 
    spread( type, cover )
  
  ggplot( avg_cover, aes( x = pred_cover, y = obs_cover )) + geom_point() + facet_grid(Group ~ Period )  + geom_abline(aes(intercept = 0, slope = 1))
  
  annual_cover <- 
    test %>%
    gather( type , cover, pred_cover, obs_cover ) %>% 
    group_by( Period, year, type )  %>% 
    summarise( cover = mean(cover) )
  
  ggplot( annual_cover, aes( x = year, y = cover, linetype = type )) + geom_line() 
  
  # ------------------------------------------------------------------------------------------------------ #
  # true_cover <- data.frame(year = gdl$year3, quad = gdl$quad3 , Period = gdl$Period3, Treatment = gdl$treat3, Group = gdl$gid3, X =  gdl$X3)
  # 
  # true_cover$size <- true_cover$X*gdl$Xscale + gdl$Xcenter # re-scale size based on mean and sd 
  # 
  # true_cover$size <- 100*(exp( as.numeric( true_cover$size) ) /10000) # convert to % cover
  # 
  # true_cover <- 
  #   true_cover %>% 
  #   group_by(Period, Treatment, year, quad ) %>% 
  #   summarise( obs_cover = sum(as.numeric( size ) ) ) %>% 
  #   group_by( Period, Treatment , year ) %>% 
  #   summarise( obs_cover = mean(obs_cover))
  # 
  # last_cover$Period <- as.numeric( factor( last_cover$Period ) )
  # 
  # trueCov <- rbind(true_cover, last_cover)

  
  
  plot_cover <-
    left_join( trueCov, cover, all.y = TRUE ) %>%
    mutate( pred_cover = ifelse(Treatment != 1 & year == 2011, obs_cover, pred_cover )) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 2007, obs_cover, pred_cover )) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 1929, obs_cover, pred_cover)) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 1945, obs_cover, pred_cover)) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 1949, obs_cover, pred_cover)) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 1954, obs_cover, pred_cover)) %>%
    gather( stat, val,  pred_cover, obs_cover) 
  
  plot_cover$Treatment <- factor(plot_cover$Treatment, labels = c('Control', 'Drought', 'Irrigation'))
  plot_cover$Period <- factor( plot_cover$Period, labels = c('Historical', 'Modern') )
  
  years <- expand.grid( Treatment = c('Control', 'Drought', 'Irrigation'), stat = c('obs_cover', 'pred_cover'), year = min(plot_cover$year):max(plot_cover$year) )
  years$Period <- 'Historical'
  years$Period[years$year > 1958 ] <- 'Not monitored'
  years$Period[years$year > 2006 ] <- 'Modern'
    
  plot_cover <- merge( years , plot_cover , all.x = T)
  
  print( 
    ggplot( subset( plot_cover, Period == "Historical"), aes( x = year, y =  val, color = Treatment, linetype = stat)) +
      geom_line() +
      scale_color_manual(values = my_colors ) + 
      scale_x_continuous()
  )
  
  print( 
    ggplot( subset( plot_cover, Period == "Modern"), aes( x = year, y =  val, color = Treatment, linetype = stat)) +
      geom_line() +
      scale_color_manual(values = my_colors ) + 
      scale_x_continuous()
    )
  
  
  dev.off()
} 



