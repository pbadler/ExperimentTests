
# call from removal_analysis_wrapper.r

if(exists("trtTests")==F)  read.csv("treatment_test_results.csv")

# plotting function
plotEffects <- function(doStage){
  ii<-which(trtTests$stage==doStage)
  myLims <- c(min(trtTests$CI.02.5[ii],na.rm=T),max(trtTests$CI.97.5[ii],na.rm=T))
  plot(c(1:4),trtTests$effect[ii],ylim=myLims,xlab="",ylab="",xaxt="n",pch=16)
  axis(1,at=c(1:4),labels=c("ARTR","HECO","POSE","PSSP"),las=2)
  abline(h=0,col="gray")
  arrows(x0=c(1:4),y0=trtTests$CI.02.5[ii],x1=c(1:4),y1=trtTests$CI.97.5[ii],
         code=3,angle=90,length=0.05) 
}

png("treatment_tests.png",height=3,width=7,res=400,units="in")

  par(mfrow=c(1,3),tcl=-0.2,mgp=c(2,0.5,0),mar=c(4,2,3,1),oma=c(2,2,0,0))
  
  plotEffects(doStage="survival")
  mtext("(A) Survival",side=3,line=1,adj=0,cex=0.9)
  
  plotEffects(doStage="growth")
  mtext("(B) Growth",side=3,line=1,adj=0,cex=0.9)
  
  #remove ARTR from recruitment plot--error bars too big
  trtTests[which(trtTests$species=="ARTR" & trtTests$stage=="recruitment"),3:5] <- NA
  plotEffects(doStage="recruitment")
  mtext("(C) Recruitment",side=3,line=1,adj=0,cex=0.9)
  
  mtext("Species",side=1,line=0,outer=T)
  mtext("Removal effect",side=2,line=0,outer=T)

dev.off()


