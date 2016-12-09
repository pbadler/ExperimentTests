

spp 
vr
lambda 
sd 


dat <- readRDS('data/temp_data/modified', vr, 'data_lists_for_stan.RDS')[[spp]]

if ( vr != 'recruitment' )  { 
  mydat <- data.frame( year = dat$year ,  C = dat$C, X = dat$X, Y=dat$Y, gid=dat$gid, yid = dat$yid  ) 
}else if( vr == 'recruitment'){ 
  mydat <- data.frame( year = dat$year ,  C = dat$C, Y=dat$Y, gid=dat$gid, yid = dat$yid, parents1 = dat$parents1, parents2 = dat$parents2  ) 
} 

year <- unique( dat$year )

for( i in 1:length(year)){ 
  
  year_out <- year[i]
  mydat_out <- mydat[mydat$year == year_out, ]
  mydat_in  <- mydat[mydat$year != year_out, ]

  mydat_in$C <- scale(mydat_in$C)
  mydat_out$C <- scale(mydat_out$C, attr(mydat_in$C, 'scaled:center', 'scaled:scale'))

  mydat_in <- as.list( mydat_in)
  mydat_in$N <- length(mydat$Y)
  mydat_in$nyrs <- length(mydat$yid)
  mydat_in$Nspp <- ncol(mydat$) 
  mydat_in$

  
  mydat_in <- as.list(mydat_in) 
  mydat_out <- as.list(mydat_out)

  names(mydat_in)
  names(mydat_out)
  
  # make new data list 
  # run stan 
  
}

  
  
  
  
