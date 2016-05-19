
# call from removal_analysis_wrapper.r

distWts <- read.csv(paste0(root,"/driversdata/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))
#distWts <- read.csv(paste0(root,"/driversdata/data/idaho/speciesdata/IdahoDistanceWeights.csv"))

# kernel figure
png("CompKernels.png",height=3,width=4.5,units="in",res=400)
par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,1,1))
matplot(distWts$midRings[1:12],distWts[1:12,1:4],xlab="Distance (cm)",ylab="Weight",
    type="l",lty=1,lwd=1.5,col=c("black","darkgray","dodgerblue3","red3"))
legend("topright",c("ARTR","HECO","POSE","PSSP"),lwd=1.5,col=c("black","darkgray","dodgerblue3","red3"),
       bty="n")
dev.off()

# make distWts constant over first annulus (0 to 2 cm)
distWts <- rbind(distWts,c(1,1,1,1,0))
distWts <- distWts[order(distWts$midRings),]

Wfuns <- list(4) # store approximate spline functions
zc <- numeric(4) # store distance at which spline function goes to zero
W.constant <- rep(0,4) # store constant for heterospecific mean field approximation
for(i in 1:4){
  
  Wfuns[[i]] <- approxfun(x=distWts$midRings, y=distWts[,i], rule=2) # W as a continuous function

  tmp <- which(distWts[,i]==0)
  if(length(tmp)>0){
    zc[i] <- distWts$midRings[min(tmp)]
  }else{
    zc[i] <- 150  # max distance
  }
  
  # if W is continuous
  W.constant[i] <- 2*pi*integrate(f=function(z) z*Wfuns[[i]](z), lower=0, upper=zc[i])$value
  
}
names(Wfuns) <- names(distWts[,1:4])

# take a look at the functions
# x<-seq(0,25,0.1)
# y<-Wfuns[[1]](x)
# plot(x,y,type="l")
# for(i in 2:4){
#   lines(x,Wfuns[[i]](x))
# }

