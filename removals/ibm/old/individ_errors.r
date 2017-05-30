# analyze individual level errors

errors1=plants
errors2=plants[,c("quad","species","year","trackID","logarea")]
names(errors2)[5] = c("logarea.obs")
errors2$year=errors2$year-1
errors1=merge(errors1,errors2,all.x=T)

pdf("obsVpred_by_plant.pdf",height=8,width=8)

par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
for(i in 1:4){
  tmp=subset(errors1, species==sppList[i])
  plot(tmp$logarea.pred~tmp$logarea.obs,main=sppList[i])
  abline(0,1,lty="dashed")
  #abline(lm(tmp$logarea.pred~tmp$logarea.obs))

}
mtext("Observed log size",1,outer=T,line=0.5)
mtext("Predicted log size",2,outer=T,line=0.5)

dev.off()


