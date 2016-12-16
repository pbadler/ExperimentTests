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

generate_cover_predictions <- function( spp, model ) { 
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[spp]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[spp]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[spp]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_', model, '_fit.RDS'), full.names = TRUE)
  
  s <- rstan::extract(readRDS(m[3]), c('muhat2'))$muhat2
  g <- rstan::extract(readRDS(m[1]), c('muhat3'))$muhat3
  r <- rstan::extract(readRDS(m[2]), c('lambda_pred2'))$lambda_pred2
  
  df <- read.csv(paste0( 'data/temp_data/', spp, '_growth_and_survival_cleaned_dataframe.csv'))
  rdf <- read.csv(paste0( 'data/temp_data/', spp, '_recruitment_cleaned_dataframe.csv'))
  
  if (model == 'treatment'){  #### Treatment model fit is very bad so just sample with replacement from real values
    r <- r 
    r[] <- sample(rdf$Y, length(r), replace = T)
  }
  
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
    summarize( total_size = sum( size ))                # sum area of all plants per quad 
  
  rdf <- 
    rdf %>% 
    gather( simulation, rec_area, matches('^X[0-9]')) 
  
  predicted_cover <- merge(rdf[ , c('year', 'Treatment', 'simulation', 'rec_area', 'quad')], cover1, all.x = T)
  
  predicted_cover$total_size <- ifelse( is.na( predicted_cover$total_size), 0, predicted_cover$total_size )
  
  predicted_cover$area <- predicted_cover$rec_area + predicted_cover$total_size
  
  predicted_cover$model <- model
  
  saveRDS( predicted_cover, paste0( 'output/ibm/simulations/', spp, '_one_step_ahead_', model, '_model_quadrat_cover.RDS'))
  
  return(predicted_cover)
}


#
load('analysis/figure_scripts/my_plotting_theme.Rdata') 

years <- expand.grid(year = 1925:2016, Treatment = c('Control', 'Drought', 'Irrigation'), type = c('observed', 'predicted'))
years$Period[ years$year > 2010 ] <- 'Modern'
years$Period[ years$year <= 2011 ] <- 'Historical'

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
model_list <- c('climate', 'treatment', 'year_effects')

ylims <- list( c(0,25), c(0,4), c(0,4), c(0,4))
names(ylims) <- species_list
iter <- expand.grid( species = species_list, model = model_list )

for( i in 1:nrow(iter) ) { 
  
  spp <- as.character( iter$species[i])
  model <- as.character( iter$model[i] )
  
  pred_cover <- generate_cover_predictions(spp, model)
  
  # get list of all quads observed each year ------------------------------------------------------------# 
  all_quad_years1 <- pred_cover %>% dplyr::select(Treatment, quad, year ) %>% distinct() # first year of transition
  all_quad_years2 <- all_quad_years1
  all_quad_years2$year <- all_quad_years2$year + 1 # also need cover predictions for each quad in the second year 
  all_quad_years <- unique( rbind( all_quad_years1, all_quad_years2 )) # both years 
  
  # get observed cover --------------------------------------------------------------------------------- #   
  df <- read.csv(paste0('data/temp_data/', spp, '_growth_and_survival_cleaned_dataframe.csv'))
  df <- df[, c('year', 'Treatment', 'quad', 'obs_id', 'logarea.t0')]
  df$area <- exp( df$logarea.t0 )
  qarea <- df %>% group_by( year, Treatment, quad ) %>% summarise( area = sum( area ))
  qarea <- merge( all_quad_years1, qarea , all.x = T) 
  qarea$area <- ifelse(is.na(qarea$area), 0, qarea$area )  # assign quadrats without plants zero cover 
  
  # get second year of observed cover, not in survival data -----------------------------------------------#  
  df <- read.csv(paste0('data/temp_data/', spp, '_growth_and_survival_cleaned_dataframe.csv'))
  df <- df[, c('year', 'Treatment', 'quad', 'obs_id', 'logarea.t1')]
  df$year <- df$year + 1  # using next years cover so add one to year 
  df$area <- exp( df$logarea.t1 ) # get cover of surviving plants in next year 
  qarea2 <- df %>% filter( !is.na(area)) %>% group_by( year, Treatment, quad ) %>% summarise( area = sum(area)) 

  # get recruit area ---------------------------------------------------------------------------------------# 
  df <- read.csv(paste0('~/driversdata/data/idaho/speciesData/', spp, '/', spp, '_genet_xy.csv'))
  df$quad <- as.numeric( str_extract(df$quad, '[0-9]+'))
  df$year <- df$year + 1900  
  df <- df %>% filter( age == 1 ) %>% group_by( quad, year ) %>% summarise( area = sum(area ) ) # get area of recruits 
  
  df2 <- read.csv(paste0( '~/driversdata/data/idaho_modern/speciesData/', spp, '/', spp, '_genet_xy.csv'))
  df2$quad <- as.numeric( str_extract(df2$quad, '[0-9]+'))
  df2 <- df2 %>% filter( age == 1 ) %>% group_by( quad, year ) %>% summarise( area = sum(area ) ) # get area of recruits 
  df2 <- subset(df2, year == 2016) # only need 2016 
  
  df <- rbind( df, df2) # last year recruits 
  
  # merge recruits with second year of cover   
  qarea2 <- merge( qarea2, df , by  = c('year', 'quad' ), all.x = T ) 
  qarea2$area.x [ is.na(qarea2$area.x) ] <- 0 
  qarea2$area.y [ is.na( qarea2$area.y ) ] <- 0 
  qarea2$area <- qarea2$area.x + qarea2$area.y # sum up area 
  
  qarea2 <-  merge(all_quad_years2, qarea2, by = c('Treatment', 'year', 'quad'), all.x = T)
  qarea2$area <- ifelse( is.na(qarea2$area), 0 , qarea2$area )
  
  # merge second year with first year 
  qarea <- merge( qarea, qarea2[, c('year', 'quad', 'Treatment', 'area')], by = c('quad', 'year', 'Treatment' ), all.x = T, all.y = T)
  qarea$area <- ifelse( is.na( qarea$area.x ), qarea$area.y, qarea$area.x)  # fill in missing area for second year 
  
  write.csv(qarea, paste0( 'output/ibm/simulations/', spp, '_observed_cover_per_quadrat.csv'), row.names = F)
  
  # merge observed and predicted cover ---------------------------------------------------
  pred_cover$year <- pred_cover$year + 1  # predicted cover is for the following year 
  
  plot_df <- merge( qarea[ , c('area', 'quad', 'year', 'Treatment') ], pred_cover[, c('area', 'quad', 'year', 'Treatment', 'simulation')], by = c('year','quad', 'Treatment') , all.x = T, all.y = T)
  
  # merge_with quad_info to get groups 
  quad_info <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')
  quad_info$quad <- as.numeric(str_extract(quad_info$quad, '[0-9]+'))
  plot_df <- merge( plot_df, quad_info)
  
  predicted_cover <- 
    plot_df %>% 
    filter( Group != 'P2') %>%
    group_by( simulation, year, Treatment  ) %>% 
    filter( !is.na( area.y) ) %>% 
    summarise( cover = 100*mean(area.y)/10000 ) %>% 
    group_by( year, Treatment ) %>% 
    summarise(predicted = mean(cover), ucl = quantile( cover, 0.75), lcl = quantile( cover, 0.25) )
  
  observed_cover <- 
    plot_df %>% 
    group_by( year, quad, Treatment ) %>% 
    summarise( area.x = unique(area.x) ) %>% 
    group_by( year, Treatment) %>% 
    summarise( observed = 100*mean(area.x)/10000 ) 
  
  plot_df <- 
    merge( observed_cover, predicted_cover, all.x = T) %>% 
    mutate( predicted = ifelse(is.na(predicted), observed, predicted))
  
  plot_df$Period <- ifelse( plot_df$year < 2011 , 'Historical', 'Modern')
  
  plot_df <- 
    plot_df %>% 
    gather( type, cover, observed, predicted ) %>% 
    mutate ( ucl = ifelse( is.na(ucl), cover, ucl )) %>% 
    mutate( lcl = ifelse (is.na(lcl), cover, lcl ))

  plot_df <- merge( years, plot_df, all.x = T, all.y = T)
  
  plot_df <- 
    plot_df %>% 
    filter( !(Treatment != 'Control' & Period == 'Historical')) %>% 
    mutate( Period = ifelse(is.na(Period), 'Historical', Period )) %>%
    mutate( Treatment = ifelse(Period == 'Historical', 'Control' , as.character(Treatment))) %>% 
    mutate( Treatment2 = ifelse( Period == 'Historical', 'Training Data', Treatment )) 
  
  plot_df$Treatment <- factor(plot_df$Treatment2, levels = c('Training Data', 'Control', 'Drought', 'Irrigation'), ordered = T)
  
  p1 <- 
    ggplot( plot_df, aes( x = year, y =  cover, fill = Treatment, color = Treatment, linetype = type, shape = type, ymax = ucl, ymin = lcl)) +
      geom_line()  +
      #geom_ribbon(aes(ymax = ucl, ymin = lcl, color = NA), alpha = 0.05) + 
      scale_color_manual(values = my_colors) + 
      scale_fill_manual( values = my_colors) +
      ylab( 'cover (%)' ) +  
      my_theme + 
      scale_y_continuous(limits = ylims[[spp]]) + 
      guides( linetype=guide_legend(title=NULL), color = guide_legend(title = NULL)) 
  
  pdf( paste0( 'figures/predictions/', spp  , '_', model, '_model_predicted_cover2.pdf' ), height = 8, width = 8) 
  
  print( 
    p1 +       
      scale_x_continuous(breaks = seq(1925,1960,5), limits = c(1928, 1959)) + 
      ggtitle(paste('Historical cover of', spp )) 
  )
  
  # fix break in line 

  plot_df <-
    rbind( plot_df  %>% 
             filter( !(year == 2011 & Period == 'Historical' )) , 
      plot_df %>% 
        filter( year == 2011 & Treatment == 'Control') %>% 
        mutate( Treatment = 'Training Data')
    ) 

  print(
    p1 %+% plot_df + 
      #geom_vline(aes(xintercept= 2011), alpha = 0.5, linetype = 4) + 
      scale_x_continuous(breaks = c(2007:2016), limits = c(2007, 2016)) + 
      ggtitle(paste('Modern cover of', spp ))
  )
  
  dev.off()
} 

