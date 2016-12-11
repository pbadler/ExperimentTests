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
  return(out )
}


#
setwd('~/Documents/ExperimentTests/precip/')

my_colors <- c('#1b9e77', '#d95f02', '#7570b3')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

ylims <- list( c(0,40), c(0,7.5), c(0,7.5), c(0,7.5))
i = 1

for( i in 1:4) {  
  spp   <- species_list[i]  
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[i]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[i]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[i]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_climate_fit.RDS'), full.names = TRUE)
  
  s1 <- rstan::extract(readRDS(m[3]), 'mu')$mu
  s2 <- rstan::extract(readRDS(m[3]), 'muhat')$muhat
  s <- cbind(s1, s2)
  
  g <- rstan::extract(readRDS(m[1]), 'muhat3')$muhat3
  
  g[ exp(g)/10000 > 1 ]  <-  log(10000)  # assign plants larger than the entire plot to the plot size 
  
  r1 <- rstan::extract(readRDS(m[2]), 'lambda')$lambda
  r2 <- rstan::extract(readRDS(m[2]), 'lambda_pred')$lambda_pred
  r <- cbind( r1, r2 )

  apply(r1, 2, quantile)
  apply(r2, 2, quantile)
  max(r1)
  min(r1)
  
  apply(r2, 2, median)
  colMeans(r2)
  
  rsize <- simulate_recruit_size( r, spp = spp)
  
  length( rsize [ exp( rsize ) /10000 > 1 ]  
  
  nsims <- dim(s)[1]
  nobs  <- dim(s)[2]
  
  # combine survival, growth and recruitment
  g <- (100*exp( g ))/10000 
  size <- s*exp(g)
  
  rsize <- (100*rsize)/10000
  
  size <- data.frame( year = gdl$year3, Treatment = gdl$treat3 , quad = gdl$quad3, group = gdl$gid3 , t(size))
  rsize  <- data.frame( year = rdl$year2, Treatment = rdl$treat2 , quad = rdl$quad2, group = rdl$gid2, t(rsize))
  
  cover1 <-  
    size %>% 
    gather( iteration, size , starts_with ('X')) %>% 
    group_by( year, quad, iteration ) %>% 
    summarise( size_cover = sum( size ))
  
  cover2 <- 
    rsize %>% 
    gather( iteration, recruit_cover, starts_with('X'))
  
  cover <- merge( cover2, cover1, all.x = T)
  
  cover$size_cover [ is.na( cover$size_cover) ] <- 0 # assign size cover to zero when no plants were present 
  
  cover$total_cover <- cover$size_cover + cover$recruit_cover 
  
  hist( cover$total_cover )
  
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
    mutate( Treatment = as.numeric(factor(Treatment))) %>%
    mutate( year = ifelse(year < 100, year + 1900, year)) %>%
    group_by(Period, Treatment, year, quad, Group ) %>%
    summarise( totCover = (100*totCover)/10000) %>% 
    group_by( Treatment, Period ) %>% 
    filter( year == max(year)) %>% 
    mutate ( quad = as.numeric(str_extract(quad, '[0-9]+') ))
  
  last_cover$Group <- as.numeric( last_cover$Group )
  
  # ------------------------------------------------------------------------------------------------------ #
  true_cover <- data.frame(year = gdl$year3, quad = gdl$quad3 , Treatment = gdl$treat3, Group = gdl$gid3, X =  gdl$X3)
  true_cover$size <- scale(true_cover$X, gdl$Xcenter, gdl$Xscale)
  true_cover$size <- exp( as.numeric( true_cover$size) )
  
  true_cover <- 
    true_cover %>% 
    group_by( year, Treatment, quad, Group ) %>% 
    summarise( totCover = sum(as.numeric( size ) ) )

  true_cover$Period <- NA
  true_cover$Period[ true_cover$year < 2006  ] <- 'Historical'
  true_cover$Period[ true_cover$year > 2006  ] <- 'Modern'
  
  hist(true_cover$totCover)
  hist(last_cover$totCover)
  
  true_cover$totCover <- (true_cover$totCover*10000)/100 
  
  trueCov <- rbind(true_cover, last_cover)
  
  
  cover <- 
    cover %>% 
    group_by( year, quad, group ) %>% 
    summarise( mcover = mean(total_cover), lci75 = quantile( total_cover, 0.125), uci75 = quantile( total_cover, 0.875))
  
  %>% 
    mutate( mcover = (mcover*100) / 10000 , lci75 = 100*lci75/10000, uci75 = uci75*100/10000)
  
  
  merge( trueCov, cover )
  
  plot_cover[[i]] <-
    merge( pred_cover[[i]], trueCov , all.y = TRUE) %>%
    mutate( pred_cover = ifelse(Treatment != 1 & year == 2011, true_cov, pred_cover )) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 2007, true_cov, pred_cover )) %>%
    gather( stat, val,  pred_cover, true_cov) 
  
  plot_cover[[i]] <- merge( years, plot_cover[[i]], all.x = TRUE)

  plot_cover[[i]]$Treatment <- factor(plot_cover[[i]]$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  pdf( paste0( 'figures/predictions/', spp  , 'lppd_predicted_cover.pdf' ), height = 8, width = 8) 
  
  print( 
    ggplot( subset( plot_cover[[i]], Period == "Historical"), aes( x = year, y =  val, color = Treatment, linetype = stat)) +
      geom_line() +
      scale_color_manual(values = my_colors ) + 
      ylim( ylims[[i]]) + 
      scale_x_continuous()
  )
  
  print( 
    ggplot( subset( plot_cover[[i]], Period == "Modern"), aes( x = year, y =  val, color = Treatment, linetype = stat)) +
      geom_line() +
      scale_color_manual(values = my_colors ) + 
      ylim( ylims[[i]] ) + 
      scale_x_continuous()
  )
  
  dev.off()
} 



