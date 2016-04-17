
# call from removal_analysis_wrapper.r

distWts <- read.csv(paste0(root,"/driversdata/data/idaho/speciesdata/IdahoDistanceWeights.csv"))

Wfuns <- list(4) # store approximate spline functions
zc <- numeric(4) # store distance at which spline function goes to zero
W.constant <- numeric(4) # store constant for heterospecific mean field approximation
for(i in 1:4){
  Wfuns[[i]] <- approxfun(x=distWts$midRings, y=distWts[,i], rule=2) 
  
  tmp <- which(distWts[,i]==0)
  if(length(tmp)>0){
    zc[i] <- distWts$midRings[min(tmp)]
  }else{
    zc[i] <- 150  # max distance
  }
  
  W.constant[i] <- 2*pi*integrate(f=function(z) z*Wfuns[[i]](z), lower=0, upper=zc[i])$value
  
}
names(Wfuns) <- names(distWts[,1:4])



