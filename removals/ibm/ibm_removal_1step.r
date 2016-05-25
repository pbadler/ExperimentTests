
# Predict survival and growth of each plant and recruitment in each quadrat*year
# PBA 5/25/2016

outfile="ibm/simulations1step/ObsPred_1step.csv"

# FORMAT PARAMETERS ------------------------------------------------
Nspp=length(sppList)

curDir <- getwd()
Nyrs <- 30
# set up survival parameters and function
source("survival/import2ibm_1step.r")
# set up growth parameters and function
source("growth/import2ibm_1step.r")
# set up recruitment parameters and function
source("recruitment/import2ibm_1step.r")
setwd(curDir)

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
  
  # drop rainfall treatment plots
  ii <- which(D2$Treatment=="Drought" | D2$Treatment=="Irrigation")
  D2 <- D2[-ii,]

  # account for removal in baseline years and take out removed plants
  if(doSpp!="ARTR"){
    ii <- which(D2$year>=2011 & D2$Treatment=="No_shrub")
    D2$W.ARTR[ii] <- 0
    ii <- which(D2$Treatment=="No_grass")
    D2<-D2[-ii,]
  }else{
    ii <- which(D2$year>=2011 & D2$Treatment=="No_grass")
    D2$W.HECO[ii] <- 0 ; D2$W.POSE[ii] <- 0 ; D2$W.PSSP[ii] <- 0
    ii <- which(D2$Treatment=="No_shrub")
    D2<-D2[-ii,]
  }

  D2$area <- exp(D2$logarea)
  D2$doSpp <- doSpp
  D2 <- subset(D2,year>2010)
  
  plants <- rbind(plants,D2)
  
}

# get quadrat group and treatment
quad.info <- unique(plants[,c("quad","Group","Treatment")],MARGIN=2)
# merge in quad codes
tmp<-data.frame(Group=c("E1","P1","P10E1","P1E1","P2","P7E1"),GroupCode=c(1:6))
quad.info <- merge(quad.info,tmp)
# add GroupCode to plants data frame too
plants <- merge(plants,tmp)

# aggregate to quadrat and year
cov.obs <- aggregate(plants$area,by=list(species=plants$doSpp,quad=plants$quad,year=plants$year),FUN=sum)
names(cov.obs)[4] <- "obs"
cov.obs$obs <- cov.obs$obs/100 # convert to % cover
cov.obs <- reshape(cov.obs,idvar=c("quad","year"),direction="wide",timevar="species")
cov.obs[is.na(cov.obs)] <- 0
cov.obs <- cov.obs[,c(1,2,6,3,4,5)] # reorder columns

# get 2015 quadrat cover totals (these are not in the survival data file)
tmp <- read.csv("QuadYearCover.csv")
tmp <- subset(tmp, year==2015)
# set removed spp to zero
tmp$cover[tmp$species=="Artemisia tripartita" & tmp$Treatment=="No_shrub"] <- 0
tmp$cover[tmp$species!="Artemisia tripartita" & tmp$Treatment=="No_grass"] <- 0
tmp <- tmp[,c("quad","year","species","cover")] # drop Treatment and Group columns
tmp <- reshape(tmp, idvar=c("quad","year"),timevar="species",direction="wide")
names(tmp) <- names(cov.obs)
cov.obs <-rbind(cov.obs,tmp)

# GET PREDICTIONS -------------------------------------------------------

# SURVIVAL AND GROWTH
plants$surv.prob <- plants$surv.prob.trt <- NA
plants$logarea.pred <- plants$logarea.pred.trt <- NA
W.index <- grep("W.",names(plants))

for(k in 1:dim(plants)[1]){
  
  doYr <- which(Spars$yrList==plants$year[k])
  doSpp <- which(sppList==plants$species[k])
  
  #ignore treatment effects
  plants$surv.prob[k]=survive(Spars,doSpp=doSpp,doGroup=plants$GroupCode[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment="Control")
  plants$logarea.pred[k] <- grow(Gpars,doSpp=doSpp,doGroup=plants$GroupCode[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment="Control")
  
  #use treatment effects when appropriate
  plants$surv.prob.trt[k]=survive(Spars,doSpp=doSpp,doGroup=plants$GroupCode[k],
      doYear=doYr,sizes=plants$logarea[k],crowding=plants[k,W.index],Treatment=plants$Treatment[k])
  plants$logarea.pred.trt[k] <- grow(Gpars,doSpp=doSpp,doGroup=plants$GroupCode[k],
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

# loop through observed cover matrix
out.recruit <- out.recruit.trt <- data.frame(quad=cov.obs$quad,year=as.numeric(cov.obs$year),matrix(NA,dim(cov.obs)[1],4))
names(out.recruit)[3:6] <- sppList
names(out.recruit.trt)[3:6] <- sppList
out.recruit <- subset(out.recruit,year<2015); out.recruit.trt <- subset(out.recruit.trt,year<2015)
for(k in 1:dim(out.recruit)[1]){
    
    totArea <- cov.obs[which(cov.obs$year==out.recruit$year[k] & cov.obs$quad==out.recruit$quad[k]),3:6]
    doYr<- which(Spars$yrList==cov.obs$year[k])
    qI <- which(quad.info$quad==cov.obs$quad[k])
    doGroup <- quad.info$GroupCode[qI]
    Trt <- quad.info$Treatment[qI]
    
    # no treatment effects
    out.recruit[k,3:6]=recruit(Rpars,totArea=as.numeric(cov.obs[k,3:6]),
          doGroup=doGroup,doYear=doYr,Treatment="Control")
    
    # with trtEffects as appropriate
    out.recruit.trt[k,3:6]=recruit(Rpars,totArea=as.numeric(cov.obs[k,3:6]),
          doGroup=doGroup,doYear=doYr,Treatment=Trt)
    
} # next k

# ADD RECRUIT COVER TO SURVIVAL*GROWTH COVER

#reorder columns
tmp <- sort(names(cov.pred)[3:NCOL(cov.pred)],index.return=T)$ix
cov.pred <- cov.pred[,c(1,2,(2+tmp))]

# make sure rows are ordered
cov.pred <- cov.pred[order(cov.pred$quad,cov.pred$year),]
out.recruit <- out.recruit[order(out.recruit$quad,out.recruit$year),]
out.recruit.trt <- out.recruit.trt[order(out.recruit.trt$quad,out.recruit.trt$year),]
cbind(cov.pred[,c(1,2)],out.recruit[,c(1,2)],out.recruit.trt[,c(1,2)]) # check rows

# now add together
cov.pred[,3:6] <- cov.pred[,3:6] + out.recruit[,c(3:6)]
cov.pred[,7:10] <- cov.pred[,7:10] + out.recruit.trt[,c(3:6)]

# FORMAT OUTPUT -------------------------------------------------------

# merge observed and predicted cover
cov.pred$year <- cov.pred$year + 1
output <- merge(cov.obs,cov.pred,all=T)

# add in quad.info
output <- merge(output,quad.info)

write.table(output,outfile,row.names=F,sep=",")

rm(plants)
