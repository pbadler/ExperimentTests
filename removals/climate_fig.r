
dat <- read.csv("../data/idaho_modern/climateData/monthly_climate_from_daily.csv")

dat$climYr <- ifelse(dat$month>9,dat$year+1,dat$year)

dat <- subset(dat,climYr>1926) # remove incomplete year

annualPpt <- aggregate(dat[,c("TPPT")],by=list(dat$climYr),FUN="sum")
names(annualPpt) <- c("climYr","precip")
annualPpt[,2] <- annualPpt[,2]*25.4 # convert in to mm
meanPpt <- mean(annualPpt[,2],na.rm=T)
ppt.90 <- quantile(annualPpt[,2],c(0.05,0.95))
modernPpt <- subset(annualPpt,climYr>2010)

annualTemp <- aggregate(dat[,c("TMEAN")],by=list(dat$climYr),FUN="mean")
names(annualTemp) <- c("climYr","temperature")
annualTemp[,2] <- (annualTemp[,2]-32)*5/9 # convert F to C
meanTemp <- mean(annualTemp[,2],na.rm=T)
temp.90 <- quantile(annualTemp[,2],c(0.05,0.95),na.rm=T)
modernTemp <- subset(annualTemp,climYr>2010)

# make figure
png("climate.png",height=3.5,width=8,res=400,units="in")

par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

plot(modernPpt,type="o",pch=16,ylim=c(150,500),xlab="Year",ylab="Annual precipitation (mm)")
abline(h=meanPpt,col="dodgerblue3",lwd=1)
abline(h=ppt.90[1],col="dodgerblue3",lwd=1,lty="dashed")
abline(h=ppt.90[2],col="dodgerblue3",lwd=1,lty="dashed")
mtext("(a)",side=3,line=0.5,adj=0)

plot(modernTemp,type="o",pch=16,ylim=c(4,9.5),xlab="Year",ylab="Mean temperature (C)")
abline(h=meanTemp,col="red3",lwd=1)
abline(h=temp.90[1],col="red3",lwd=1,lty="dashed")
abline(h=temp.90[2],col="red3",lwd=1,lty="dashed")
mtext("(b)",side=3,line=0.5,adj=0)

dev.off()
