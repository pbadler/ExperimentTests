
# call from removal_analysis_wrapper.r

sppList=c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
dataDir <- paste(root,"/driversdata/data/idaho_modern/",sep="")


# import data and calculate treatment trends ######################################

ptsD<-read.csv(paste(dataDir,"allrecords_density.csv",sep=""))
ptsD$density <- 1
trts<-read.csv(paste(dataDir,"quad_info.csv",sep=""))

# use this to make sure we don't miss zeros
allquadyrs<-unique(ptsD[,c("quad","year")],MARGIN=2)
allquadyrs<-merge(allquadyrs,trts[,c("Group","quad","Treatment")])

# only do removal treatments
keep <- which(is.element(allquadyrs$Treatment,c("Control","No_shrub","No_grass")))
allquadyrs <- allquadyrs[keep,]

# aggregate to quadrat level
keep<-which(!is.element(ptsD$Species,sppList))
sppD<-ptsD[keep,]
sppD<-merge(sppD,allquadyrs,all.y=T)
sppD$density[is.na(sppD$density)]<-0
sppD.q<-aggregate(sppD$density,by=list(Group=sppD$Group,Treatment=sppD$Treatment,
                  quad=sppD$quad,year=sppD$year),FUN=sum)
names(sppD.q)[NCOL(sppD.q)]<-"density"


# focus on all control plots or just those in the big exclosure?
 #sppD.q <- subset(sppD.q,Group=="E1")

#calculate treatment means by year
spp.mean <- aggregate(sppD.q$density,by=list(Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=mean)
names(spp.mean)[NCOL(spp.mean)] <- "density"
spp.mean <- reshape(spp.mean,direction="wide",timevar="Treatment",idvar=c("year"))
spp.mean <- subset(spp.mean,year>2010)


# figures ########################################################################

myCols<-c("black","darkgoldenrod","darkgreen")

#1. Average cover treatment and year
matplot(spp.mean$year,spp.mean[,2:4],xlab="Year",ylab="Cover (%)",
        type="o",pch=16,lty=1,col=myCols)
