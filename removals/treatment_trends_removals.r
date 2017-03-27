
# call from removal_analysis_wrapper.r

sppList=c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
dataDir <- paste(root,"/ExperimentTests/data/idaho_modern/",sep="")


# import data and calculate treatment trends ######################################

covD<-read.csv(paste(dataDir,"allrecords_cover.csv",sep=""))
trts<-read.csv(paste(dataDir,"quad_info.csv",sep=""))

# use this to make sure we don't miss zeros
allquadyrs<-unique(covD[,c("quad","year")],MARGIN=2)
tmp<-expand.grid(species=sppList,year=sort(unique(covD$year)))
allquadyrs<-merge(allquadyrs,tmp,all=T)
allquadyrs<-merge(allquadyrs,trts[,c("Group","quad","Treatment")])

# only do removal treatments
keep <- which(is.element(allquadyrs$Treatment,c("Control","No_shrub","No_grass")))
allquadyrs <- allquadyrs[keep,]

# aggregate to quadrat level
keep<-which(is.element(covD$species,sppList))
sppD<-covD[keep,]
sppD<-merge(sppD,allquadyrs,all.y=T)
sppD$area[is.na(sppD$area)]<-0
sppD$area<-sppD$area*100
sppD.q<-aggregate(sppD$area,by=list(species=sppD$species,Group=sppD$Group,Treatment=sppD$Treatment,
                  quad=sppD$quad,year=sppD$year),FUN=sum)
names(sppD.q)[NCOL(sppD.q)]<-"cover"
write.csv(sppD.q,"QuadYearCover.csv",row.names = FALSE)  # save these 

# focus on all control plots or just those in the big exclosure?
# sppD.q <- subset(sppD.q,Group=="E1")

#calculate treatment means by year
spp.mean <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=mean)
names(spp.mean)[NCOL(spp.mean)] <- "cover"
spp.mean <- reshape(spp.mean,direction="wide",timevar="Treatment",idvar=c("species","year"))
spp.mean <- subset(spp.mean,year>2010)
spp.mean <- spp.mean[order(spp.mean$species,spp.mean$year),]

#calculate standard errors by year
spp.sd <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=sd)
spp.N <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=length)
spp.se <- spp.sd
names(spp.se)[NCOL(spp.se)] <- "StErr"
spp.se$StErr <- spp.se$StErr/sqrt(spp.N$x)
rm(spp.sd);rm(spp.N)
spp.se <- reshape(spp.se,direction="wide",timevar="Treatment",idvar=c("species","year"))
spp.se <- subset(spp.se,year>2010)
spp.se <- spp.se[order(spp.se$species,spp.se$year),]

# calculate year-to-year log changes
tmp <- sppD.q
tmp$year <- tmp$year + 1
names(tmp)[which(names(tmp)=="cover")] <- "lag.cover"
logChange <- merge(sppD.q,tmp)
logChange$pcgr <- log(logChange$cover/logChange$lag.cover)
logChange$pcgr[logChange$pcgr==Inf] <- NA
logChange$pcgr[logChange$pcgr==-Inf] <- NA
mean.change <- aggregate(logChange$pcgr,by=list(species=logChange$species,Treatment=logChange$Treatment,
                  year=logChange$year),FUN=mean,na.rm=T)
names(mean.change)[NCOL(mean.change)] <- "pcgr"
mean.change <- subset(mean.change,year>2010)
mean.change <- reshape(mean.change,direction="wide",timevar="Treatment",idvar=c("species","year"))
mean.change <- mean.change[order(mean.change$species,mean.change$year),]

# calculate deviations from pretreatment year
sppD.q.2011<-subset(sppD.q,year==2011) #get pre-treatment year
names(sppD.q.2011)[NCOL(sppD.q.2011)]<-"cover.2011"
i<-which(names(sppD.q.2011)=="year")
sppD.q.2011<-sppD.q.2011[,-i]
sppD.q<-subset(sppD.q,year>2010)
sppD.q<-merge(sppD.q,sppD.q.2011)
sppD.q$coverDiff<-sppD.q$cover-sppD.q$cover.2011
spp.mean.diff<-aggregate(sppD.q$coverDiff,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=mean)
names(spp.mean.diff)[NCOL(spp.mean.diff)]<-"coverDiff"
spp.mean.diff <- reshape(spp.mean.diff,direction="wide",timevar="Treatment",idvar=c("species","year"))
spp.mean.diff <- spp.mean.diff[order(spp.mean.diff$species,spp.mean.diff$year),]


# statistical tests ####################################################

#library(lme4)

dARTR <- subset(logChange,species=="Artemisia tripartita" & !is.na(pcgr) & Treatment!="No_shrub")
dARTR$year <- as.factor(dARTR$year)
mARTR <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dARTR)

dHECO <- subset(logChange,species=="Hesperostipa comata" & !is.na(pcgr) & Treatment!="No_grass")
dHECO$year <- as.factor(dHECO$year)
mHECO <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dHECO)

dPOSE <- subset(logChange,species=="Poa secunda" & !is.na(pcgr) & Treatment!="No_grass")
dPOSE$year <- as.factor(dPOSE$year)
mPOSE <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dPOSE)

dPSSP <- subset(logChange,species=="Pseudoroegneria spicata" & !is.na(pcgr) & Treatment!="No_grass")
dPSSP$year <- as.factor(dPSSP$year)
mPSSP <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dPSSP)

texreg(list(mARTR,mHECO,mPOSE,mPSSP), ci.force=TRUE,caption="Cover change models",
      caption.above=TRUE,file=statsOutput)


# figures ########################################################################

trtLabels<-substr(x=names(spp.mean)[3:5],start=7,stop=nchar(names(spp.mean)[3:5]))
myCols<-c("black","darkgoldenrod","darkgreen")

#1. Average cover treatment and year
png("treatment_trends_cover.png",height=6,width=8,units="in",res=400)
  
  par(mfrow=c(2,2),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0),tcl=-0.2)

  # mean cover
  for(doSpp in sppList){
    tmp.mean<-subset(spp.mean,species==doSpp)
    tmp.se<-subset(spp.se,species==doSpp)
    if(doSpp==sppList[1]){
      my.y <- c(0,36)
    }else{
      my.y<- c(0,4.2)
    }
    matplot(tmp.mean$year,tmp.mean[,3:5],ylim=my.y,type="o",xlab="",ylab="",pch=16,lty="solid",
            col=myCols,main=doSpp,font.main=4,lwd=2)
    for(iTrt in 1:3){
      arrows(x0=tmp.mean$year,y0=tmp.mean[,2+iTrt]+tmp.se[,2+iTrt],x1=tmp.mean$year,y1=tmp.mean[,2+iTrt]-tmp.se[,2+iTrt],
             col=,myCols[iTrt],length=0.02,code=3,angle=90)
    }
    if(doSpp==sppList[2])  legend("topright",c("Control","Grass removal","Shrub removal"),pch=16,lty="solid",col=myCols,bty="n")
      
  }
  mtext("Mean cover (%)",side=2,line=0.5,outer=T,cex=1.2)
  mtext("Year",side=1,line=0.5,outer=T,cex=1.2)
  
dev.off()

#1. Average cover treatment and year (this version drops the data for the removal species)
png("treatment_trends_cover_simple.png",height=2.75,width=8,units="in",res=400)
  
  par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0),tcl=-0.2)

  # mean cover
  for(doSpp in sppList){
    tmp.mean<-subset(spp.mean,species==doSpp)
    my.y <-c(0,max(tmp.mean[,3:5]))
    if(doSpp=="Artemisia tripartita"){
      tmp.mean$cover.No_shrub <- NA 
    }else{
      tmp.mean$cover.No_grass <- NA
    }
    matplot(tmp.mean$year,tmp.mean[,3:5],ylim=my.y,type="o",xlab="",ylab="",pch=16,lty="solid",
            col=myCols,main=doSpp,font.main=4,lwd=2)
    if(doSpp==sppList[1]) {
      legend("topright",c("Control","Grass removal","Shrub removal"),pch=16,lty="solid",col=myCols,bty="n")
      mtext("Mean cover (%)",side=2,line=2,outer=F)
    }
  }
  
  mtext("Year",side=1,line=1,outer=T)
  
dev.off()

#2. Log cover change by treatment and year 
png("treatment_trends_logChange.png",height=3,width=8.5,units="in",res=400)
  
  par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0),tcl=-0.2,lwd=1)

  # log change
  myLims <- c(-1,1)
  for(doSpp in sppList){
    tmp.mean<-subset(mean.change,species==doSpp )
    tmp.mean[1,3:5] <- NA  # get rid of control plot value in 2011
    #remove irrelevant removal treatments
    if(doSpp=="Artemisia tripartita"){
      tmp.mean$pcgr.No_shrub <- NA
    }else{
      tmp.mean$pcgr.No_grass <- NA
    }
    matplot(tmp.mean$year,tmp.mean[,3:5],type="o",xlab="",ylab="",pch=16,lty="solid",
            col=myCols,ylim=myLims,main=doSpp, font.main=4)
    abline(h=0,col="gray")
    
    if(doSpp==sppList[1]) {
      legend("topright",c("Control","Grass removal","Shrub removal"),pch=16,lty="solid",col=myCols,bty="n")
      mtext("Mean log change",side=2,line=2,outer=F)
    }
    
  }
  mtext("Year",side=1,line=0.25,outer=T)

dev.off()


