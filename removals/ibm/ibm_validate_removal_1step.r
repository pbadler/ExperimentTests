# Individually-based model for a multiple species,
# with density-dependence and explicit space (random
# spatial pattern)

# Simulate one quadrat through time

# PBA  8-28-09

# qName = "Q1"
Ngroups=6
startYr=2011
# doGroup=1  # E1 exclosure
L=100 # dimension of square quadrat (cm)
expand=1  # 1 = 1x1 m^2, 2 = 2x2m^2, etc
#sppList=c("ARTR","HECO","POSE","PSSP")
myCol=c("black","gold1","blue","red")
minSize=0.25
maxSize=c(8000,500,500,500)

if(trtEffects==F){
 outfile1=paste("simulations1step/",qName,"_validation_cov_removals_noTrt.csv",sep="")
 outfile2=paste("simulations1step/",qName,"_validation_den_removals_noTrt.csv",sep="")
}else{
 outfile1=paste("simulations1step/",qName,"_validation_cov_removals_Trt.csv",sep="")
 outfile2=paste("simulations1step/",qName,"_validation_den_removals_Trt.csv",sep="") 
}

# FORMAT PARAMETERS ------------------------------------------------
Nspp=length(sppList)

curDir <- getwd()
Nyrs <- 30
# set up survival parameters and function
source("survival/import2ibm_1step.r")
# set up growth parameters and function
source("growth/import2ibm_deterministic.r")
# set up recruitment parameters and function
source("recruitment/import2ibm_deterministic.r")
setwd(curDir)

# model spatial group variation (or not)
if(!is.na(doGroup)){
  Spars$intcpt=Spars$intcpt+Spars$intcpt.gr[doGroup,]
  Gpars$intcpt=Gpars$intcpt+Gpars$intcpt.gr[doGroup,]
  Rpars$intcpt.yr=Rpars$intcpt.yr+matrix(Rpars$intcpt.gr[doGroup,],Nyrs,Nspp,byrow=T)
}

# FUNCTIONS---------------------------------------------------------
library(boot)
library(mvtnorm)
library(msm)
source("survival/fetchSurvData.r")

# GET OBSERVED DATA  -------------------------------------------
Nspp=length(sppList)
plants=NULL
neighborList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)
dataDir2 <- paste(root,"/driversdata/data/idaho_modern",sep="")

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  D2 <- fetchSdat(doSpp=doSpp,speciesList=neighborList,datadir=dataDir2,distWts=dists)

  # merge in treatment data
  tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
  tmp <- tmp[,c("quad","Treatment")]
  D2 <- merge(D2,tmp, all.x=T)

  # account for removal in baseline years
  if(doSpp!="ARTR"){
    ii <- which(D2$year>=2011 & D2$Treatment=="No_shrub")
    D2$W.ARTR[ii] <- 0
  }else{
    ii <- which(D2$year>=2011 & D2$Treatment=="No_grass")
     D2$W.HECO[ii] <- 0 ; D2$W.POSE[ii] <- 0 ; D2$W.PSSP[ii] <- 0
  }
  
  # clean up dataset 
  if(doSpp=="ARTR"){
    keep <- which(is.element(D2$Treatment,c("Control","No_grass")))
  }else{
    keep <- which(is.element(D2$Treatment,c("Control","No_shrub")))
  }
  D2 <- D2[keep,]
  D2 <- subset(D2, year>2010) # start with 2011

  D2$area <- exp(D2$logarea)
  D2$doSpp <- doSpp
  
  plants <- rbind(plants,D2)
  
}

# aggregate to quadrat and year
cov.obs <- aggregate(plants$area,by=list(species=plants$doSpp,quad=plants$quad,year=plants$year),FUN=sum)
names(cov.obs)[4] <- "cover"
cov.obs$cover <- cov.obs$cover/100 # convert to % cover
cov.obs <- reshape(cov.obs,idvar=c("quad","year"),direction="wide",timevar="species")
cov.obs[is.na(cov.obs)] <- 0
cov.obs <- cov.obs[,c(1,2,6,3,4,5)] # reorder columns



# GET PREDICTIONS -------------------------------------------------------

# survival and growth
plants$surv.prob <- NA; plants$logarea.pred <- NA
W.index <- grep("W.",names(plants))[1:4]  # TODO: import allcov and allpts coefficients

for(k in 1:dim(plants)[1]){
  
  doYr <- which(Spars$yrList==plants$year[k])
  doSpp <- which(sppList==plants$species[k])
  plants$surv.prob[k]=survive(Spars,doSpp=plants$species[k],doGroup=plants$Group[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=matrix(plants[k,W.index],1,4))
}

# arrays to store results
N=matrix(0,(simYrs+1),Nspp)
N[1,1:4]=as.numeric(obsN[1,2:5])
A=matrix(0,(simYrs+1),Nspp)
A[1,1:4]=as.numeric(obsA[1,2:5])

for(tt in 1:simYrs){
  
  if(is.null(init.plants[[tt]])){ # are there any plants?
    
    A[tt+1,] <- 0
    N[tt+1,] <- 0
    
  }else{
    
    # initialize with N.init plants of size.init for each species
    plants=init.plants[[tt]]
    lastID=max(plants[,5])
    
    # draw year effects
    doYr=doYrList[tt]  # no year effects
    
    nextplants=plants
    
    # recruitment
    newplants=recruit(Rpars,sizes=plants[,2],spp=plants[,1],doGroup=doGroup,doYear=doYr,lastID=lastID,L,expand)
    
    W=getCrowding(plants,L,expand)
    
    for(ss in 1:Nspp){
      if(N[tt,ss]>0){ # make sure spp ss is not extinct
        
        # growth
        newsizes=grow(Gpars,doSpp=ss,doGroup=doGroup,doYear=doYr,sizes=plants[,2],crowding=W)
        if(sum(newsizes==Inf)>0) browser()
        if(is.na(sum(newsizes))) browser()
        
        # survival, uses same W as growth            
        live=survive(Spars,doSpp=ss,doGroup=doGroup,doYear=doYr,sizes=plants[,2],crowding=W)
        
        # put it all together
        tmp=which(plants[,1]==ss)  # only alter plants of focal spp        
        nextplants[tmp,2]=newsizes[tmp]*live[tmp]   #update with G and S
        
      } # end if no plants
    } # next ss  
    
    nextplants=nextplants[nextplants[,2]>0,]    # remove dead plants 
    nextplants=rbind(nextplants,newplants)     # add recruits
    
    if(dim(nextplants)[1]==0) break()  # end simulation
    
    # output cover and density
    tmp=aggregate(nextplants[,2],by=list(nextplants[,1]),FUN=sum)
    A[tt+1,tmp[,1]]=tmp[,2]
    tmp=aggregate(rep(1,dim(nextplants)[1]),by=list(nextplants[,1]),FUN=sum)
    N[tt+1,tmp[,1]]=tmp[,2]
    
    # since we are re-initializing every year, no need for this
    #      plants=nextplants
    #      lastID=max(plants[,5])  
    
  } # end is.null(plants) loop
  
  print(tt);flush.console() 
  
} # next tt 

print(paste(qName," complete",sep=""))
flush.console()

# format output

year=c(calYrList[1],1+calYrList)

A=data.frame(cbind(year,A))
names(A)[2:dim(A)[2]]=paste(sppList,"pred",sep="")
output1=merge(obsA,A,all.x=T)

N=data.frame(cbind(year,N))
names(N)[2:dim(N)[2]]=paste(sppList,"pred",sep="")
output2=merge(obsN,N,all.x=T)

# par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0))
# matplot(output1[,1],output1[,2:NCOL(output1)],type="o",
#   col=myCol,lty=c(rep("solid",Nspp),c(rep("dashed",Nspp))),
#   pch=c(rep(16,Nspp),rep(1,Nspp)),xlab="Year",ylab="Cover")
# matplot(output2[,1],output2[,2:NCOL(output2)],type="o",
#   col=myCol,lty=c(rep("solid",Nspp),c(rep("dashed",Nspp))),
#   pch=c(rep(16,Nspp),rep(1,Nspp)),xlab="Year",ylab="Density")

write.table(output1,paste0("ibm/",outfile1),row.names=F,sep=",")
write.table(output2,paste0("ibm/",outfile2),row.names=F,sep=",")

