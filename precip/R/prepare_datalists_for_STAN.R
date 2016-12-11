######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

rm(list = ls() )

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

detachAllPackages()

source('R/climate/ExtractData_3Runs.R')
source('R/climate/aggregate_spot_VWC.R')
source('R/climate/soilMoistureTreatmentEffects.R')
source('R/climate/aggregate_VWC_data.R')
source('R/get_all_demographic_data.R')
source('R/climate/make_climate_variables.R')
source('R/climate/prepare_climate_covariates.R')

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(stringr)

make_data_list <- function( x, vr ) { 
  
  x$X <- scale(x$X)
  Xcenter <- attr(x$X, 'scaled:center')
  Xscale  <- attr(x$X, 'scaled:scale')
  x$X <- as.numeric(x$X)
  
  W <- x[ , grep ( '^W\\.', names( x))]
  W <- as.matrix( W )[,1:4] # big four competition effects
  W <- scale(W)
  x$W <- W
  Wcenter <- attr( W, 'scaled:center')
  Wscale  <- attr( W, 'scaled:scale' )

  x$C <- as.matrix( x[ ,  grep( '^C\\.', names(x) ) ] )
  colnames(x$C) <- str_replace( colnames( x$C ) , '^C\\.' , '')
  
  x$C <- scale( x$C )
  
  Ccenter <- attr(x$C, 'scaled:center')
  Cscale  <- attr(x$C, 'scaled:scale')
  # x$Chold <- scale( x$Chold, x$Ccenter, x$Cscale)
  # x$C2    <- scale( x$C2, x$Ccenter, x$Cscale)
  
  x$treat <- as.numeric(factor( x$Treatment))
  x$tm <- model.matrix.lm( ~ x$Treatment )[, -1 ]  # drop intercept 
  x$tm <- cbind ( x$tm , x$tm*x$X ) 
  colnames(x$tm) <- c('Drought' , 'Irrigation', 'Droughtxlogarea.t0', 'Irrigationxlogarea.t0')
  
  x$gm <- model.matrix.lm( ~ x$Group ) 
  x$spp <- as.numeric(factor( x$species))
  
  saveRDS(x, file =  paste0( 'data/temp_data/', unique( x$species) ,'_', vr, '_cleaned_dataframe.RDS'))
  
  x <- 
    x %>% 
    dplyr::select( year, treat, Period, spp, quad, trackID, W, obs_id, C, yid, gid, X, Y, gm, tm ) %>% 
    rename( trackid = trackID)

  mylist <- split( x , x$Period)
  
  mylist <- lapply( mylist, as.list )
  
  mylist$all <- as.list( x )
  
  mylist <- lapply( mylist, function( y ) { y$nyrs = nlevels(factor(y$yid)); y } )
  mylist <- lapply( mylist, function( y ) { y$nT = ncol(y$tm); y } )
  mylist <- lapply( mylist, function( y ) { y$G = nlevels(factor(y$gid)); y } )
  mylist <- lapply( mylist, function( y ) { y$N = length(y$Y); y})
  
  names(mylist$Modern) <- paste( names(mylist$Modern) ,'hold', sep = '' ) 
  names(mylist$all) <- paste( names(mylist$all), 2, sep = '')
  
  mylist <- unlist(mylist, recursive = F, use.names = T)
  
  names( mylist ) <- str_replace( names(mylist) , '^.*\\.', '') # clean up names 
  
  mylist$Period <- as.numeric( factor( mylist$Period ) )
  mylist$Periodhold <- as.numeric( factor(mylist$Periodhold))
  mylist$Period2 <- as.numeric( factor(mylist$Period2))
  
  mylist$Ccenter <- Ccenter
  mylist$Cscale <- Cscale
  
  mylist$Wcenter <- Wcenter
  mylist$Wscale <- Wscale 

  mylist$Xcenter <- Xcenter 
  mylist$Xscale <- Xscale   
  
  mylist$Wcovs <- ncol ( mylist$W )
  mylist$Covs <- ncol ( mylist$C ) 
  mylist$spp <- unique(x$spp)
  
  mylist$tau_beta <- 10 
  
  return( mylist ) 
  
} 


add_survival_data <- function( growth , survival ) { 

  survival$W2 <- scale( (survival$W2*survival$Wscale + survival$Wcenter), growth$Wcenter, growth$Wscale ) 
  survival$C2 <- scale( (survival$C2*survival$Cscale + survival$Ccenter), growth$Ccenter, growth$Cscale ) 
  survival$X2 <- scale( (survival$X2*survival$Xscale + survival$Xcenter), growth$Xcenter, growth$Xscale ) 
  
  survival$X2 <- as.numeric( survival$X2 )
  
  cover_list <- survival [ grep ( '*.2$', names( survival ) ) ] 
  names(cover_list )  <- str_replace(names(cover_list) , '2$', '3')
  
  out <- c( growth, cover_list)
  
  return( out ) 

} 

make_data_list_recruitment <- function( x, vr ) { 
  
  x$obs_id <- as.numeric(row.names(x))
  x$yid <- as.numeric(factor(x$year))
  x$quad <- as.numeric(factor(x$quad))
  
  x$gid <- as.numeric(x$Group)
  x$treat <- as.numeric(factor( x$Treatment))
  x$gm <- model.matrix.lm(~ x$Group)
  x$tm <- model.matrix.lm(~ x$Treatment)[ , -1]
  
  x$C <- as.matrix( x[ ,  grep( '^C\\.', names(x) ) ] )
  colnames(x$C) <- str_replace( colnames( x$C ) , '^C\\.' , '')
  
  x$spp <- as.numeric(factor( x$species))
  
  x$parents1 <- as.matrix( x[ , grep( '^cov\\.', names(x))] )
  x$parents2 <- as.matrix( x[ , grep( '^Gcov\\.', names(x))] )
  
  saveRDS(x, paste0( 'data/temp_data/', unique(x$species ) ,'_', vr, '_cleaned_dataframe.RDS'))
  
  x <- 
    x %>% 
    dplyr::select( Y, year, treat, Period, spp, quad, obs_id, C, yid, gid, parents1, parents2, gm, tm ) 
  
  mylist <- split( x , x$Period)
  
  mylist <- lapply( mylist, as.list )
  
  mylist$all <- as.list( x )
  
  mylist <- lapply( mylist, function( y ) { y$nyrs = nlevels(factor(y$yid)); y } )
  mylist <- lapply( mylist, function( y ) { y$nT = ncol(y$tm); y } )
  mylist <- lapply( mylist, function( y ) { y$G = nlevels(factor(y$gid)); y } )
  mylist <- lapply( mylist, function( y ) { y$N = length(y$Y); y})
  
  names(mylist$Modern) <- paste( names(mylist$Modern) ,'hold', sep = '' ) 
  names(mylist$all) <- paste( names(mylist$all), 2, sep = '')
  
  mylist <- unlist(mylist, recursive = F, use.names = T)

  names( mylist ) <- str_replace( names(mylist) , '^.*\\.', '') # clean up names 
  
  mylist$Period <- as.numeric( factor( mylist$Period ) )
  mylist$Periodhold <- as.numeric( factor(mylist$Periodhold))
  mylist$Period2 <- as.numeric( factor(mylist$Period2))
  
  mylist$C <- scale( mylist$C )
  mylist$Ccenter <- attr(mylist$C, 'scaled:center')
  mylist$Cscale  <- attr(mylist$C, 'scaled:scale')
  mylist$Chold <- scale( mylist$Chold, mylist$Ccenter, mylist$Cscale)
  mylist$C2    <- scale( mylist$C2, mylist$Ccenter, mylist$Cscale)
  
  mylist$Nspp <- ncol ( mylist$parents1 )
  mylist$Covs <- ncol ( mylist$C ) 
  mylist$spp <- unique( x$spp ) 
  mylist$tau_beta <- 10 
  
  return( mylist ) 
} 

# -- select covariates -------------------------------------------------------------------#

clim_vars <- c('VWC.sp.l', 
               'VWC.sp.0', 
               'VWC.sp.1',
               'VWC.su.l', 
               'VWC.su.0', 
               'VWC.su.1',
               'VWC.f.l', 
               'VWC.f.0', 
               'VWC.f.1',
               'T.sp.1', 
               'T.sp.0',
               'T.sp.l',
               'T.su.1', 
               'T.su.0', 
               'T.su.l',
               'T.f.1', 
               'T.f.0', 
               'T.f.l', 
               'T.w.1',
               'T.w.0', 
               'T.w.l')                     

clim_file <- 'all_clim_covs.RDS'
data_path <- 'data/temp_data'

species <- c('ARTR', 'HECO', 'POSE', 'PSSP')
spp_num <- as.factor(species )

clim <- readRDS(paste0( data_path, '/', clim_file ) )  

clim <- clim[ ,c('Treatment', 'year', clim_vars)]
names ( clim ) [ grep( '^(P\\.)|(VWC\\.)|(T\\.)',  names( clim ) ) ] <- paste0 ( 'C.', names( clim ) [ grep( '^(P\\.)|(VWC\\.)|(T\\.)', names(clim ) ) ] )

out <- list()
i = 1
for(i in 1:length( species )) { 
  spp <- species[i]
  
  gdat <- readRDS( paste0(data_path, '/', spp, '_growth.RDS'))
  sdat <- readRDS( paste0(data_path, '/', spp, '_survival.RDS'))
  rdat <- readRDS( paste0( data_path, '/', spp, '_recruitment.RDS'))
  
  sdat$logarea.t0 <- sdat$logarea
  sdat <- sdat [ sdat$year %in% gdat$year, ] 
  
  df <- merge( sdat, gdat, all.x = T)
  
  rm(gdat, sdat)
  
  df <- df %>% arrange( year, Group, quad, trackID )
  df$obs_id <- as.numeric( row.names(df) )
  
  clim <- clim[ complete.cases(clim), ] 
  
  df <- merge(df, clim)
  rdf <- merge( rdat, clim ) 
  
  df <- 
    df %>% 
    mutate( yid = as.numeric( factor(year)), 
            gid = as.numeric( factor(Group)), 
            quad = as.numeric(factor(quad)) ) 
  
  # split into data lists 
  
  df$X <- df$logarea.t0
  
  growth <- df[ !is.na(df$logarea.t1), ]
  survival <- df[ !is.na(df$survives), ]
  
  survival$Y <- survival$survives
  growth$Y <- growth$logarea.t1
  
  growth <- make_data_list( x = growth, 'growth' )  
  survival <- make_data_list(x = survival, 'survival' )
  recruitment <- make_data_list_recruitment( x = rdf, 'recruitment' ) 
  
  growth <- add_survival_data(growth, survival )

  out[[i]] <-  list( growth, recruitment, survival ) 
  
  rm(df, growth, survival, recruitment)
}

ARTR <- out[[1]][[1]]
HECO <- out[[2]][[1]]
POSE <- out[[3]][[1]]
PSSP <- out[[4]][[1]]

saveRDS( list ( ARTR = ARTR, HECO = HECO, POSE = POSE, PSSP = PSSP ), 'data/temp_data/growth_data_lists_for_stan.RDS')

ARTR <- out[[1]][[2]]
HECO <- out[[2]][[2]]
POSE <- out[[3]][[2]]
PSSP <- out[[4]][[2]]

saveRDS( list ( ARTR = ARTR, HECO = HECO, POSE = POSE, PSSP = PSSP ), 'data/temp_data/recruitment_data_lists_for_stan.RDS')

ARTR <- out[[1]][[3]]
HECO <- out[[2]][[3]]
POSE <- out[[3]][[3]]
PSSP <- out[[4]][[3]]

saveRDS( list ( ARTR = ARTR, HECO = HECO, POSE = POSE, PSSP = PSSP ), 'data/temp_data/survival_data_lists_for_stan.RDS')
