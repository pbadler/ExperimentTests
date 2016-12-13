rm(list = ls())

var_table <- read.csv('output/selected_climate_covariates.csv')

for( i in 1:nrow(var_table)){ 
  
  spp    <- var_table$species[i]
  vr <- var_table$vital_rate[i]
  covars     <- strsplit( as.character( var_table$covars[i] ) , ',')[[1]]
  
  dl <- readRDS(paste0('data/temp_data/', vr, '_data_lists_for_stan.RDS'))
  
  dat         <- dl[[spp]]
  
  # make interactions --------------------------------------------- # 
  if( vr != 'recruitment') { 
    ifx <- dat$C*as.numeric( dat$X )
    colnames(ifx) <- paste(colnames(dat$C), 'logarea.t0', sep = 'x')
    C  <- cbind( dat$C, ifx) 
    dat$C <- as.matrix(C)

    ifx <- dat$C2*as.numeric(dat$X2)
    colnames(ifx) <- paste(colnames(dat$C2), 'logarea.t0', sep = 'x')
    dat$C2  <- cbind( dat$C2, ifx)
    
    ifx <- dat$Chold*as.numeric( dat$Xhold )
    colnames(ifx) <- paste(colnames(dat$Chold), 'logarea.t0', sep = 'x')
    dat$Chold  <- cbind( dat$Chold, ifx)
    
    if ( vr == 'growth' ) { 
      ifx <- dat$C3*as.numeric(dat$X3)
      colnames(ifx) <- paste(colnames(dat$C3), 'logarea.t0', sep = 'x')
      dat$C3  <- cbind( dat$C3, ifx ) 
      
      dat$C3 <- dat$C3[, covars]
      dat$Covs3 <- ncol( dat$C3)
    }
  }
   
  dat$C       <- dat$C[ , covars, drop = F] 
  dat$C2      <- dat$C2[, covars, drop = F]
  dat$Chold   <- dat$Chold[, covars, drop = F]
    
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



