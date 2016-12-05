rm(list = ls())

library(rstan)

var_table <- read.csv('output/selected_climate_covariates.csv')
i = 1
for( i in 1:nrow(var_table)){ 
  
  spp    <- var_table$species[i]
  vr <- var_table$vital_rate[i]
  covars     <- strsplit( as.character( var_table$covars[i] ) , ',')[[1]]
  
  dl <- readRDS(paste0('data/temp_data/', vr, '_data_lists_for_stan.RDS'))
  
  dat         <- dl[[spp]]
  
  # make interactions --------------------------------------------- # 
  if( vr != 'recruitment') { 
    ifx <- dat$C*dat$X
    colnames(ifx) <- paste(colnames(dat$C), 'logarea.t0', sep = 'x')
    C  <- scale( cbind( dat$C, ifx) ) 
    dat$C <- C
    dat$Ccenter <- attr(C, 'scaled:center')
    dat$Cscale  <- attr(C, 'scaled:scale')
    
    ifx <- dat$C2*dat$X2
    colnames(ifx) <- paste(colnames(dat$C2), 'logarea.t0', sep = 'x')
    dat$C2  <- scale( cbind( dat$C2, ifx), dat$Ccenter, dat$Cscale )
    
    ifx <- dat$Chold*dat$Xhold
    colnames(ifx) <- paste(colnames(dat$Chold), 'logarea.t0', sep = 'x')
    dat$Chold  <- scale( cbind( dat$Chold, ifx), dat$Ccenter, dat$Cscale)
    
    if ( vr == 'growth' ) { 
      ifx <- dat$C3*dat$X3
      colnames(ifx) <- paste(colnames(dat$C3), 'logarea.t0', sep = 'x')
      dat$C3  <- scale( cbind( dat$C3, ifx ) , dat$Ccenter, dat$Cscale)
      
      dat$C3 <- dat$C3[, covars]
      dat$Covs3 <- ncol( dat$C3)
    }
  }
   
  dat$C       <- dat$C[ , covars] 
  dat$C2      <- dat$C2[, covars]
  dat$Chold   <- dat$Chold[, covars]
  
  dat$Ccenter <- dat$Ccenter[covars]
  dat$Cscale  <- dat$Cscale[covars]
  
  dat$Ccenter <- dat$Ccenter[ !is.na(dat$Ccenter)]
  dat$Cscale  <- dat$Cscale[ !is.na(dat$Cscale)]
    
  dat$Covs    <- ncol(dat$C)
  
  if(file.exists(paste0('data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))) { 
    mod_dl <- readRDS(paste0('data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))
    mod_dl[[spp]] <- dat
  }else{ 
    mod_dl <- dl 
    mod_dl[[spp]] <- dat
  }
  
  saveRDS(mod_dl, paste0( 'data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))
}



