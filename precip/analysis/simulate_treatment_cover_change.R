rm(list =ls())
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)

# functions ---------------------------------------------------------------------------- 

get_survival_size_covariates <- function( sdl, gdl ) {
 
  return( 
    list(  
      id_table = data.frame( year = sdl$year2, Treatment = sdl$treat2, yid = sdl$yid2, quad = sdl$quad2, trackID = sdl$trackid2),
      yid      = sdl$yid2,
      nyrs     = sdl$nyrs2,
      X        = sdl$X2, 
      C        = as.matrix(sdl$C2),
      W        = sdl$W2,
      gm       = sdl$gm2,
      N        = sdl$N2,
      Xcenter  = gdl$Xcenter,  # use growth Yscale and Ycenter 
      Xscale   = gdl$Xscale 
      )
  )
} 

get_growth_size_covariates <- function( sdl, gdl ) {
  
  return( 
    list(  
      id_table = data.frame( year = gdl$year3, Treatment = gdl$treat3, yid = gdl$yid3, quad = gdl$quad3, trackID = gdl$trackid3),
      yid      = gdl$yid3,
      nyrs     = gdl$nyrs3,
      X        = gdl$X3, 
      C        = as.matrix(gdl$C3),
      W        = gdl$W3,
      gm       = gdl$gm3,
      N        = gdl$N3,
      Xcenter  = gdl$Xcenter,  # use growth Yscale and Ycenter 
      Xscale   = gdl$Xscale 
    )
  )
} 



get_recruitment_covariates <- function( rdl ) {
  
  return( 
    list( 
      id_table = data.frame( year = rdl$year2, Treatment = rdl$treat2, yid = rdl$yid2, quad = rdl$quad2), 
      yid      = rdl$yid2, 
      nyrs     = rdl$nyrs2,
      spp_id   = rdl$spp,
      C        = as.matrix(rdl$C2),
      parents1 = rdl$parents12,
      parents2 = rdl$parents22,
      gm       = rdl$gm2,
      N        = rdl$N2
    )
  )
} 


simulate_survival <- function( pars , test_dat ){ 
  
  nyrs  <- test_dat$nyrs
  N <- test_dat$N
  C <- test_dat$C
  W <- test_dat$W
  X <- test_dat$X
  gm <- test_dat$gm
  yid <- test_dat$yid
  
  attach(pars )
  Y <- matrix( NA, nsims, N)
  
  for( j in 1:nsims){ 
    b2 <- as.matrix(b2)
    climEff  <- C%*%b2[j, ]
    gint     <- gm%*%bg[j, ]
    coverEff <- W%*%w[j, ]
    
    # draw random year effects 
    a  <- rnorm(nyrs, 0, sd = sig_a[j])
    b1 <- rnorm(nyrs, b1_mu[j], sd = sig_b1[j])
    
    mu <- coverEff
    
    for(n in 1:N){
      mu[n] <- inv.logit(gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + coverEff[n] + climEff[n])
    }
    
    #Y[j, ] <- rbinom(N, 1, mu)
    Y[j, ] <- mu
  }
  rm(pars)
  detach(pars )   
  return(Y)
}


simulate_growth <- function( pars , test_dat ){ 
  
  nyrs  <- test_dat$nyrs
  N <- test_dat$N
  C <- test_dat$C
  W <- test_dat$W
  X <- test_dat$X
  gm <- test_dat$gm
  yid <- test_dat$yid
  Xcenter <- test_dat$Xcenter
  Xscale   <- test_dat$Xscale
  
  attach(pars )
  Y <- matrix( NA, nsims, N)
  
  for( j in 1:nsims){ 
    climEff  <- C%*%b2[j, ]
    gint     <- gm%*%bg[j, ]
    coverEff <- W%*%w[j, ]
    
    # draw random year effects 
    a  <- rnorm(nyrs, 0, sd = sig_a[j])
    b1 <- rnorm(nyrs, b1_mu[j], sd = sig_b1[j])
    
    mu <- coverEff

    for(n in 1:N){
      mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + coverEff[n] + climEff[n]
    }
    
    #Y[j, ] <- rnorm(N, mu, sigma)
    Y[j, ] <- mu
  }
  
  Y <- exp( Y) # re-scale 
  
  rm(pars)
  detach(pars )   
  return(Y)
}


simulate_recruitment <- function( pars , test_dat ){ 
  
  nyrs  <- test_dat$nyrs
  N <- test_dat$N
  C <- test_dat$C
  parents1 <- test_dat$parents1
  parents2 <- test_dat$parents2
  gm <- test_dat$gm
  yid <- test_dat$yid
  spp_id <- test_dat$spp_id
  
  attach(pars )
  
  Y <- matrix( NA, nsims, N)
  
  for( j in 1:nsims){ 
    
    trueP1 <- parents1*u[j] + parents2*(1-u[j])
    trueP2 <- trueP1
    
    trueP2 <- sqrt(trueP1)
    
    climEff  <- C%*%b2[j, ]
    gint     <- gm%*%bg[j, ]
    coverEff <- trueP2%*%w[j, ]
    
    # draw random year effects 
    a  <- rnorm(nyrs, 0, sd = pars$sig_a[j])
    
    mu <- coverEff
    lambda <- coverEff
    
    for(n in 1:N){
      mu[n] <- exp(gint[n] + a[yid[n]] + coverEff[n] + climEff[n]);
      lambda[n] <- trueP1[n, spp_id]*mu[n];  
    }
    
    #Y[j, ] <- rnbinom(N, mu = lambda, size = theta[j] )
    Y[j, ] <- lambda
  }
  rm(pars)
  detach(pars )  
  return(Y)
}

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

setwd('~/Documents/ExperimentTests/precip/')

my_colors <- c('#1b9e77', '#d95f02', '#7570b3')

years <- expand.grid(year = 1925:2017, Treatment = c(1:3), stat = c('true_cov', 'pred_cover'))
years$Period[ years$year > 2006 ] <- 'Modern'
years$Period[ years$year <= 1960 ] <- 'Historical'

# output lists 
survives <- list(NA)
size     <- list(NA)
recruits <- list(NA)
rec_area <- list(NA)
k        <- list(NA)
cover    <- list(NA)
pred_cover <- list(NA)
plot_cover <- list(NA)

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

ylims <- list( c(0,40), c(0,7.5), c(0,7.5), c(0,7.5))
i = 1

for( i in 1:4) {  
  spp   <- species_list[i]  
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[i]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[i]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[i]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_climate_fit.RDS'), full.names = TRUE)
  
  spars <- rstan::extract(readRDS(m[3]), c('sig_a', 'b1_mu', 'sig_b1', 'w', 'b2', 'bg'))
  
  gpars <- rstan::extract(readRDS(m[1]), c('sig_a', 'b1_mu', 'sig_b1', 'w', 'b2', 'bg', 'sigma'))
  
  rpars <- rstan::extract(readRDS(m[2]), c('sig_a', 'w', 'b2', 'bg', 'theta', 'u'))
  
  nsims <- length(spars$sig_a)
  rpars$nsims <- gpars$nsims <- spars$nsims <- nsims
  
  sdl <-  get_survival_size_covariates(sdl, gdl)
  gdl <-  get_growth_size_covariates(sdl, gdl )
  rdl <-  get_recruitment_covariates(rdl)

  survives[[i]] <- simulate_survival( spars, sdl )
  size[[i]]     <- simulate_growth(gpars, gdl )
  recruits[[i]] <- simulate_recruitment(rpars, rdl)
  rec_area[[i]] <- simulate_recruit_size(recruits[[i]], spp )
    
  k[[i]]        <- survives[[i]]*size[[i]]  # survival by size 
  
  survives[[i]] <- data.frame(sdl$id_table , t(survives[[i]]))
  size[[i]]     <- data.frame(sdl$id_table , t(size[[i]]))
  recruits[[i]] <- data.frame(rdl$id_table , t(recruits[[i]]))
  rec_area[[i]] <- data.frame(rdl$id_table , t(rec_area[[i]]))
  k[[i]]        <- data.frame(sdl$id_table , t(k[[i]]))      
  
  k[[i]] <-   
    k[[i]] %>% 
    gather( simulation, size , starts_with('X')) %>% 
    group_by( year, Treatment, simulation, quad ) %>% 
    summarize( total_size = sum( size )) %>% 
    group_by( year, Treatment, simulation ) %>% 
    summarise( total_size = mean(total_size ))
  
  rec_area[[i]] <- 
    rec_area[[i]] %>% 
    gather( simulation, rec_area, starts_with('X')) %>% 
    group_by( year, Treatment, simulation, quad ) %>% 
    summarize( total_rec_area = sum(rec_area )) %>% 
    group_by( year, Treatment, simulation) %>% 
    summarise( total_rec_area = mean(total_rec_area ))
  
  cover[[i]] <- 
    merge( k[[i]], rec_area[[i]], all.y = TRUE) %>%  # use all recruitment years  
    mutate( total_size = ifelse(is.na(total_size), 0, total_size))  %>%  # years without observed plants get 0
    mutate( cover = 100*((total_size + total_rec_area)/10000) )
  
  pred_cover[[i]] <- 
    cover[[i]] %>% 
    mutate( year = year + 1 ) %>% # cover predictions the "Y's" are for the next year 
    group_by(Treatment, year ) %>% 
    summarise( pred_cover = median(cover), 
               ucl = quantile( cover, 0.95), 
               lcl = quantile(cover , 0.05))
  
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
    group_by(Period, Treatment, year) %>%
    summarise( true_cov = 100*(mean(totCover )/10000)) %>% 
    group_by( Treatment, Period ) %>% 
    filter( year == max(year))
  
  # ------------------------------------------------------------------------------------------------------ #
  
  sdf <- data.frame(sdl$id_table, X  = sdl$X, Xcenter =  sdl$Xcenter, Xscale = sdl$Xscale)
  sdf$Period[sdf$year < 2000 ] <- 'Historical'
  sdf$Period[sdf$year >= 2000 ] <- 'Modern'
  
  trueCov <- 
    sdf %>% 
    group_by( Period, year, Treatment, quad ) %>% 
    summarise( totCover  = sum( exp( X*Xscale + Xcenter))) %>% 
    group_by( Period, year, Treatment) %>% 
    summarise( true_cov  = 100*(mean(totCover )/10000))
  
  rm(sdl, rdl)
  
  trueCov <- rbind(trueCov, last_cover)
  
  plot_cover[[i]] <-
    merge( pred_cover[[i]], trueCov , all.y = TRUE) %>%
    mutate( pred_cover = ifelse(Treatment != 1 & year == 2011, true_cov, pred_cover )) %>%
    mutate( pred_cover = ifelse(Treatment == 1 & year == 2007, true_cov, pred_cover )) %>%
    gather( stat, val,  pred_cover, true_cov) 
  
  plot_cover[[i]] <- merge( years, plot_cover[[i]], all.x = TRUE)

  plot_cover[[i]]$Treatment <- factor(plot_cover[[i]]$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  pdf( paste0( 'figures/predictions/', spp  , '_predicted_cover.pdf' ), height = 8, width = 8) 
  
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

saveRDS(size , 'output/predicted_size.RDS')
saveRDS(survives, 'output/predicted_survival.RDS')
