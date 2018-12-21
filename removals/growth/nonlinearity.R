plot(density(sqrt(allD$W.POSE)),xlab="sqrt(W)",main="Focal species POSE",lwd=2)
lines(density(sqrt(allD$W.PSSP)),col="red",lwd=2)
lines(density(sqrt(allD$W.HECO)),col="blue",lwd=2)
lines(density(sqrt(allD$W.ARTR)),col="gold",lwd=2)
legend("topright",c("ARTR","HECO","POSE","PSSP"),
       lwd=2,col=c("gold","blue","black","red"))

m1<-lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.PSSP+ W.allcov + W.allpts +
             (logarea.t0|year),data=allD)
print(quantile(allD$W.POSE))
m1<-lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
             (logarea.t0|year),data=allD,
          subset=allD$W.PSSP < quantile(allD$W.PSSP,0.7))
m1<-lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.allcov + W.allpts +
             (logarea.t0|year),data=allD,
          subset=allD$W.PSSP < quantile(allD$W.PSSP,0.7))
