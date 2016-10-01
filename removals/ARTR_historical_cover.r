rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

sppList=c("Artemisia tripartita")
dataDir <- paste(root,"/driversdata/data/idaho/",sep="")

# import data and calculate treatment trends ######################################

covD<-read.csv(paste(dataDir,"allrecords_cover.csv",sep=""))
covD <- subset(covD,year<73)

# use this to make sure we don't miss zeros
allquadyrs<-unique(covD[,c("quad","year")],MARGIN=2)
tmp<-expand.grid(species=sppList,year=sort(unique(covD$year)))
allquadyrs<-merge(allquadyrs,tmp,all=T)

# aggregate to quadrat level
keep<-which(is.element(covD$species,sppList))
sppD<-covD[keep,]
sppD<-merge(sppD,allquadyrs,all.y=T)
sppD$area[is.na(sppD$area)]<-0
sppD$area<-sppD$area*100
sppD.q<-aggregate(sppD$area,by=list(species=sppD$species,quad=sppD$quad,year=sppD$year),FUN=sum)
names(sppD.q)[NCOL(sppD.q)]<-"cover"

#calculate treatment means by year
spp.mean <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,
                  year=sppD.q$year),FUN=mean)
names(spp.mean)[NCOL(spp.mean)] <- "cover"

#calculate treatment medians by year
spp.median <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,
                  year=sppD.q$year),FUN=median)
names(spp.median)[NCOL(spp.median)] <- "cover"

# plot all points and annual means and medians
plot(sppD.q$year, sppD.q$cover)

lines(spp.mean$year,spp.mean$cover,col="blue")
abline(h=mean(spp.mean$cover[spp.mean$year>30]),lty="dashed",col="blue")
print(mean(spp.mean$cover[spp.mean$year>30]))

lines(spp.median$year,spp.median$cover,col="red")
abline(h=mean(spp.median$cover[spp.median$year>30]),lty="dashed",col="red")
print(mean(spp.median$cover[spp.median$year>30]))
