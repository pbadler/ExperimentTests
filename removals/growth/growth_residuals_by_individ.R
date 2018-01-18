
# call this script from removal_analysis_wrapper.r

if(exists("growth_residuals")==F) stop("You need to fit growth models to get residuals")

spp_names <- c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
spp_list <- c("ARTR","HECO","POSE","PSSP")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")

# get pretreatment W's and merge to individual residuals
for(i in 1:length(spp_list)){
  doSpp <- spp_list[i]
  D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)
  if(i==1){
    growth_residuals[[i]] <- subset(growth_residuals[[1]],Treatment=="No_grass")
    D2<-subset(D2,year==2011)
    D2$W.removal <- rowSums(D2[,c("W.HECO","W.POSE","W.PSSP")])
    growth_residuals[[i]] <- merge(growth_residuals[[i]],D2[,c("quad","trackID","W.removal")],all.x=T)
  }else{
    growth_residuals[[i]] <- subset(growth_residuals[[i]],Treatment=="No_shrub")
    D2<-subset(D2,year==2011)
    D2$W.removal <- D2$W.ARTR
    growth_residuals[[i]] <- merge(growth_residuals[[i]],D2[,c("quad","trackID","W.removal")],all.x=T)
  }
}

# plot figure
png("growth_residuals_vs_Wremoval.png",height=5.5,width=8,units="in",res=400)
par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
for(i in 1:4){

  plot(growth_residuals[[i]]$W.removal,growth_residuals[[i]]$resids,xlab="",ylab="")
  abline(h=0,lty="dashed")
  mtext(side=3,spp_names[i],line=0.5,adj=0,font=3)
  
  # fit simple linear regression, display results if significant
  tmp <- lm(growth_residuals[[i]]$resids~growth_residuals[[i]]$W.removal)
  #print(summary(tmp))
  if(anova(tmp)$Pr[1] < 0.05){
    abline(tmp,col="black")
    legend("topright",paste("P =", round(anova(tmp)$Pr[1],3)),bty="n")
  }
  
}

mtext(side=1, "Pre-treatment crowding from removal species (%)",outer=T,line=0.5,cex=1.1)
mtext(side=2, "Growth residuals",outer=T,line=0.5, cex=1.1)

dev.off()
