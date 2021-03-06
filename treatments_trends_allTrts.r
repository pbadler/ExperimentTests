
rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests",sep="")); # modify as needed 

sppList=c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
dataDir <- paste(root,"/driversdata/data/idaho_modern/",sep="")

covD<-read.csv(paste(dataDir,"allrecords_cover.csv",sep=""))
trts<-read.csv(paste(dataDir,"quad_info.csv",sep=""))

# use this to make sure we don't miss zeros
allquadyrs<-unique(covD[,c("quad","year")],MARGIN=2)
tmp<-expand.grid(species=sppList,year=sort(unique(covD$year)))
allquadyrs<-merge(allquadyrs,tmp,all=T)
allquadyrs<-merge(allquadyrs,trts[,c("Group","quad","Treatment")])

# focus on exclosures
allquadyrs <- subset(allquadyrs,Group=="E1")

# calculate treatment means by year
keep<-which(is.element(covD$species,sppList))
sppD<-covD[keep,]
sppD<-merge(sppD,allquadyrs,all.y=T)
sppD$area[is.na(sppD$area)]<-0
sppD$area<-sppD$area*100
sppD.q<-aggregate(sppD$area,by=list(species=sppD$species,Treatment=sppD$Treatment,
                  quad=sppD$quad,year=sppD$year),FUN=sum)
names(sppD.q)[NCOL(sppD.q)]<-"cover"
spp.mean <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=mean)
names(spp.mean)[NCOL(spp.mean)] <- "cover"
spp.mean <- reshape(spp.mean,direction="wide",timevar="Treatment",idvar=c("species","year"))
spp.mean <- subset(spp.mean,year>2010)
spp.mean <- spp.mean[order(spp.mean$species,spp.mean$year),]

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
# statistical tests
library(lme4)

dARTR <- subset(logChange,species=="Artemisia tripartita" & !is.na(pcgr) & Treatment!="No_shrub")
dARTR$year <- as.factor(dARTR$year)
mARTR <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dARTR)
anova(mARTR)

dHECO <- subset(logChange,species=="Hesperostipa comata" & !is.na(pcgr) & Treatment!="No_grass")
dHECO$year <- as.factor(dHECO$year)
mHECO <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dHECO)
anova(mHECO)

dPOSE <- subset(logChange,species=="Poa secunda" & !is.na(pcgr) & Treatment!="No_grass")
dPOSE$year <- as.factor(dPOSE$year)
mPOSE <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dPOSE)
anova(mPOSE)

dPSSP <- subset(logChange,species=="Pseudoroegneria spicata" & !is.na(pcgr) & Treatment!="No_grass")
dPSSP$year <- as.factor(dPSSP$year)
mPSSP <- lmer(pcgr ~ Treatment + (1|quad) + (1|year),data=dPSSP)
anova(mPSSP)

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


# figures ########################################################################

trtLabels<-substr(x=names(spp.mean)[3:7],start=7,stop=nchar(names(spp.mean)[3:7]))
myCols<-c("black","red3","blue3","darkgreen","green")

#1. Average cover by treatment and year
pdf("cover_by_treatment.pdf",height=3,width=8)
par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2.5,0,0),tcl=-0.2,lwd=1)
for(doSpp in sppList){
  tmp.mean<-subset(spp.mean,species==doSpp)
  matplot(tmp.mean$year,tmp.mean[,3:7],type="o",xlab="",ylab="",pch=16,lty="solid",
          col=myCols,main=doSpp)
  if(doSpp==sppList[1]) legend("topright",trtLabels,pch=16,lty="solid",col=myCols,bty="n")
}
mtext("Year",side=1,line=1,outer=T)
mtext("Mean cover (%)",side=2,line=1,outer=T)
dev.off()

#2. Average cover deviation (w.r.t. pretreatment year)
pdf("cover_deviation.pdf",height=3,width=8)
par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2.5,0,0),tcl=-0.2)
for(doSpp in sppList){
  tmp.mean<-subset(spp.mean.diff,species==doSpp)
  matplot(tmp.mean$year,tmp.mean[,3:7],type="o",xlab="",ylab="",pch=16,lty="solid",
          col=myCols,main=doSpp)

  if(doSpp==sppList[1]) legend("bottomright",trtLabels,pch=16,lty="solid",col=myCols,bty="n")
}
mtext("Year",side=1,line=1,outer=T)
mtext("Mean cover  deviation (%)",side=2,line=1,outer=T)
dev.off()


#3. log change
# hard wire ylims
myLims <- c(-1.8,1.2)
pdf("log_change.pdf",height=3,width=8)
par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2.5,0,0),tcl=-0.2)
for(doSpp in sppList){
  tmp.mean<-subset(mean.change,species==doSpp)
  #remove irrelevant removal treatments
  if(doSpp=="Artemisia tripartita"){
    tmp.mean$pcgr.No_shrub <- NA
  }else{
    tmp.mean$pcgr.No_grass <- NA
  }
  matplot(tmp.mean$year,tmp.mean[,3:7],type="o",xlab="",ylab="",pch=16,lty="solid",
          col=myCols,main=doSpp,ylim=myLims)
  
  if(doSpp==sppList[1]) legend("bottomright",trtLabels,pch=16,lty="solid",col=myCols,bty="n")
}
mtext("Year",side=1,line=1,outer=T)
mtext("Mean cover  deviation (%)",side=2,line=1,outer=T)
dev.off()
