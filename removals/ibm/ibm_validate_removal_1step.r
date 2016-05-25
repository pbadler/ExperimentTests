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
}else{
 outfile1=paste("simulations1step/",qName,"_validation_cov_removals_Trt.csv",sep="")
}

# FORMAT PARAMETERS ------------------------------------------------
Nspp=length(sppList)

curDir <- getwd()
Nyrs <- 30
# set up survival parameters and function
source("survival/import2ibm_1step.r")
# set up growth parameters and function
source("growth/import2ibm_1step.r")
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

  D2$area <- exp(D2$logarea)
  D2$doSpp <- doSpp
  
  plants <- rbind(plants,D2)
  
}

# make quad.info data frame

# subset to > 2010

# aggregate to quadrat and year
cov.obs <- aggregate(plants$area,by=list(species=plants$doSpp,quad=plants$quad,year=plants$year),FUN=sum)
names(cov.obs)[4] <- "obs"
cov.obs$obs <- cov.obs$obs/100 # convert to % cover
cov.obs <- reshape(cov.obs,idvar=c("quad","year"),direction="wide",timevar="species")
cov.obs[is.na(cov.obs)] <- 0
cov.obs <- cov.obs[,c(1,2,6,3,4,5)] # reorder columns

# get 2015 quadrat cover totals (these are not in the survival data file)


# GET PREDICTIONS -------------------------------------------------------

# SURVIVAL AND GROWTH
plants$surv.prob <- plants$surv.prob.trt <- NA
plants$logarea.pred <- plants$logarea.pred.trt <- NA
W.index <- grep("W.",names(plants))

for(k in 1:dim(plants)[1]){
  
  doYr <- which(Spars$yrList==plants$year[k])
  doSpp <- which(sppList==plants$species[k])
  
  #ignore treatment effects
  plants$surv.prob[k]=survive(Spars,doSpp=doSpp,doGroup=plants$Group[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment="Control")
  plants$logarea.pred[k] <- grow(Gpars,doSpp=doSpp,doGroup=plants$Group[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment="Control")
  
  #use treatment effects when appropriate
  plants$surv.prob.trt[k]=survive(Spars,doSpp=doSpp,doGroup=plants$Group[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment=plants$Treatment[k])
  plants$logarea.pred.trt[k] <- grow(Gpars,doSpp=doSpp,doGroup=plants$Group[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment=plants$Treatment[k])
  
}

# multiply predicted size by survival probability to get expected area
plants$area.pred <- plants$surv.prob*exp(plants$logarea.pred)
plants$area.pred.trt <- plants$surv.prob.trt*exp(plants$logarea.pred.trt)

# aggregate predicted area to quadrat level
cov.pred <- aggregate(plants[,c("area.pred","area.pred.trt")],by=list(species=plants$doSpp,quad=plants$quad,year=plants$year),FUN=sum)
names(cov.pred)[4:5] <- c("pred","pred.trt")
cov.pred[,4:5] <- cov.pred[,4:5]/100 # convert to % cover
cov.pred <- reshape(cov.pred,idvar=c("quad","year"),direction="wide",timevar="species")
cov.pred[is.na(cov.pred)] <- 0

# RECRUITMENT

# get quadrat group and treatment
quad.info <- unique(plants[,c("quad","Group","Treatment")],MARGIN=2)

# loop through observed cover matrix
out.recruit <- out.recruit.trt <- data.frame(quad=cov.obs$quad,year=as.numeric(cov.obs$year),matrix(NA,dim(cov.obs)[1],4))
names(out.recruit)[3:6] <- sppList
names(out.recruit.trt)[3:6] <- sppList
for(k in 1:dim(cov.obs)[1]){
    
    doYr<- which(Spars$yrList==cov.obs$year[k])
    qI <- which(quad.info$quad==cov.obs$quad[k])
    doGroup <- quad.info$Group[qI]
    Trt <- quad.info$Treatment[qI]
    
    # no treatment effects
    out.recruit[k,3:6]=recruit(Rpars,totArea=as.numeric(cov.obs[k,3:6]),
          doGroup=doGroup,doYear=doYr,Treatment="Control")
    
    # with trtEffects as appropriate
    out.recruit.trt[k,3:6]=recruit(Rpars,totArea=as.numeric(cov.obs[k,3:6]),
          doGroup=doGroup,doYear=doYr,Treatment=Trt)
    
} # next k

# add recruitment to survival*growth

# format output

year=c(calYrList[1],1+calYrList)

A=data.frame(cbind(year,A))
names(A)[2:dim(A)[2]]=paste(sppList,"pred",sep="")
output1=merge(obsA,A,all.x=T)

write.table(output1,paste0("ibm/",outfile1),row.names=F,sep=",")

rm(plants)
