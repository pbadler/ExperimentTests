
# call from removal_analysis_wrapper.r

distWts <- read.csv(paste0(root,"/driversdata/data/idaho/speciesdata/IdahoDistanceWeights.csv"))
Wfuns <- list(4)
for(i in 1:4){
  Wfuns[[i]] <- approxfun(x=distWts$midRings, y=distWts[,i], rule=2)
}
names(Wfuns) <- names(distWts[,1:4])
