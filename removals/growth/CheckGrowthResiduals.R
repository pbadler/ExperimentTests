resids = (allD$logarea.t1-m1$summary.fitted.values[,1])
preds = m1$summary.fitted.values[,1]

dev.new(); par(mfrow=c(2,2)); 
eC = which(allD$Treatment=="Control");
eR = which(allD$Treatment=="No_shrub");
plot(preds[eC],resids[eC],xlab="Predicted", ylab="Observed-Predicted",ylim=range(resids)); title(main="POSE");  
points(preds[eR],resids[eR],col="red") 

scatter.smooth(preds[eC],resids[eC],main="Control"); abline(0,0,col="red",lty=2); 
scatter.smooth(preds[eR],resids[eR],main="Removal"); abline(0,0,col="red",lty=2); 
scatter.smooth(preds[eR],resids[eR],main="Removal",ylim=c(-1,1)); abline(0,0,col="red",lty=2); 
abline(h=0.2,col="blue",lty=2); 