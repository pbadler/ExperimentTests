
# call this script from removal_analysis_wrapper.r

if(exists("growth_residuals")==F) stop("You need to fit growth models to get residuals")

spp_names <- c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")

# get pretreatment cover
covD <- read.csv("C:\\Repos\\ExperimentTests\\removals\\QuadYearCover.csv")

# shrub cover
shrubD <- subset(covD, species=="Artemisia tripartita" & Treatment=="No_shrub" & year==2011)
shrubD <- shrubD[,c("quad","cover")]

# grass cover
grassD <- subset(covD,species!="Artemisia tripartita" & Treatment=="No_grass" & year==2011)
grassD <- aggregate(grassD$cover,by=list(quad=grassD$quad),FUN="sum")
names(grassD) <- c("quad","cover")

#merge in cover data
for(i in 1:4){
  if(i==1){
    growth_residuals[[i]] <- merge(growth_residuals[[i]],grassD,all.x=T)
  }else{
    growth_residuals[[i]] <- merge(growth_residuals[[i]],shrubD,all.x=T)
  }
}

# plot figure
par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,2,1),oma=c(2,2,0,0))
for(i in 1:4){
  if(i==1){
    keep <- which(growth_residuals[[i]]$Treatment=="No_grass")
  }else{
    keep <- which(growth_residuals[[i]]$Treatment=="No_shrub")
  }
  plot(growth_residuals[[i]]$cover[keep],growth_residuals[[i]]$resids[keep],xlab="",ylab="")
  abline(h=0,lty="dashed")
  mtext(side=3,spp_names[i],line=0.5,adj=0,font=3)
}
mtext(side=1, "Pre-treatment cover of removal species",outer=T,cex=1.1)
mtext(side=2, "Growth residuals",outer=T,cex=1.1)
