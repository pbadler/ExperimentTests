# Multispecies, spatially implicit IPM
# This version makes it possible to assume "no overlap"
# for intraspecific competition only or intra- and interspecific competition

#============================================================================================#
# Calculate the equilibrium areas.
#============================================================================================# 

## initial population density vector
nt=v
for(i in 1:Nspp) nt[[i]][]=0  # make sure all species = 0
for(i in init.species) nt[[i]][]=0.1  # give init.species density > 0
new.nt=nt

# set up matrix to record cover
covSave = matrix(NA,tlimit,Nspp)
covSave[1,]=sumCover(v,nt,h,A)

# set up list to store size distributions
sizeSave=list(NULL)
for(i in 1:Nspp){
  sizeSave[[i]]=matrix(NA,length(v[[i]]),(tlimit))
  sizeSave[[i]][,1]=nt[[i]]/sum(nt[[i]])
}

# initial densities 
Nsave=matrix(NA,tlimit,Nspp)
Nsave[1,]=sumN(nt,h)

yrSave=rep(NA,tlimit)
for (i in 2:(tlimit)){
  
  #draw from observed year effects
  doYear=sample(c(1:Nyrs),1)
  yrSave[i]=doYear
  
  #get recruits per area
  cover=covSave[i-1,] ; N=Nsave[i-1,]
  rpa=get.rpa(Rpars,cover,doYear)
  
  #calculate size-specific crowding
  alpha=Spars$alpha 
  for(ii in 1:Nspp){ 
    
    # first do all overlap W's
    Xbar=cover*A/N       # multiply by A to get cover back in cm^2
    varX=varN(v,nt,h,Xbar,N) 
    muW = pi*Xbar*N/(A*alpha[ii,])
    muW[is.na(muW)]=0
    Wmat[[ii]]=matrix(muW,nrow=length(v[[ii]]),ncol=Nspp,byrow=T)
    
    #now do conspecific no overlap W
    Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
    Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural")
    Wmat[[ii]][,ii]=Wrij(v.r[[ii]],ii,ii)/A
    
  }
 
  for(doSpp in 1:Nspp){  
    if(cover[doSpp]>0){    
      # make kernels and project
      K.matrix=make.K.matrix(v[[doSpp]],Wmat[[doSpp]],Rpars,rpa,Gpars,Spars,doYear,doSpp)	
      new.nt[[doSpp]]=K.matrix%*%nt[[doSpp]] 
      sizeSave[[doSpp]][,i]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
    }    
  } # next species
	
  nt=new.nt 
	covSave[i,]=sumCover(v,nt,h,A)  # store the cover as cm^2/cm^2
  Nsave[i,]=sumN(nt,h)
 	if(i%%10==0) {print(i); flush.console()}
 	if(sum(is.na(nt))>0) browser()  
 	
} # next time step

## Figures ==============================================================

par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0)) 
myCol=c("black","gold","blue","red")
#cover
boxplot(as.data.frame(100*covSave[(burn.in+1):tlimit,]),ylab="Cover (%)",names=sppList,col=myCol)
abline(h=0)
#density
boxplot(as.data.frame(Nsave[(burn.in+1):tlimit,]),ylab="Density",names=sppList,col=myCol)
abline(h=0)
# average size distribution  
plot(1,1,type="n",xlim=c(log(0.15),log(max(maxSize*2))),ylim=c(0,0.1),
  xlab="Size",ylab="Frequency")
for(i in 1:Nspp){
 lines(v[[i]],rowMeans(sizeSave[[i]][,(burn.in+1):tlimit]),col=myCol[i])
}
# example time series
matplot((burn.in+1):tlimit,100*covSave[(burn.in+1):tlimit,],type="l",col=myCol,
  xlab="Time",ylab="Cover (%)")

# # plot changes in size distribution
# X11()
# par(mfrow=c(2,2),tcl=-0.2)
# for(i in 1:Nspp) {image(sizeSave[[i]][],x=v[[i]],y=1:tlimit,
#   xlab="size",ylab="time",main=sppList[i]) }

meanCover <- colMeans(covSave[(burn.in+1):tlimit,])*100
