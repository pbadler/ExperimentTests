
### plot a size distribution after running an IPM simulation
### and compare to observed

doSpp=2

plot(1,1,type="n",xlim=c(log(0.15),log(max(maxSize*2))),ylim=c(0,0.3),
  xlab="Size",ylab="Frequency")
#simulated
lines(v[[doSpp]],rowMeans(sizeSave[[doSpp]][,(burn.in+1):tlimit]),col="red")

#observed
infile <- paste0(root,"\\ExperimentTests\\data\\idaho_modern\\speciesData\\",sppList[doSpp],"\\survD.csv")
tmp=read.csv(infile)
obsSize <- table(cut(log(tmp$area),breaks=c(v[[doSpp]]-h[doSpp],6)))
obsSize <- obsSize/sum(obsSize)
lines(v[[doSpp]],obsSize,col="black")

### compare baseline and removal predictions for quadrat level cover
# first run IBM summarize_sims1step.r
ratio=numeric(4)
for(i in 1:4){
  if(i==1) {myTrt="No_grass"}else{myTrt="No_shrub"}
  removal.base=simD[simD$Treatment==myTrt,(6+i)]
  removal.trt=simD[simD$Treatment==myTrt,(10+i)]
  tmp=removal.trt/removal.base
  ratio[i]=mean(tmp,na.rm=T)
}
print(ratio)

### POSE growth models
# first run POSEgrowth.r

# play with LMER models

#baseline model
m1.lmer <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
              (logarea.t0|year),data=allD)
summary(m1.lmer)

# try treatment x size interactions
m2.lmer <- lmer(logarea.t1~logarea.t0*Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
              (logarea.t0|year),data=allD)
summary(m2.lmer) # not significant, but close

m3.lmer<- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO  + W.POSE + W.PSSP+W.POSE*Treatment+ W.allcov + W.allpts+
              (1|Group)+(logarea.t0|year),data=allD) 
summary(m3.lmer) # SIGNIFICANT

# try fitting with weights
allD$weights <- (allD$logarea.t0+2)/(sum(allD$logarea.t0+2))
m1.weights <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
              (logarea.t0|year),weights=weights,data=allD)
summary(m1.weights) # makes treatment effect even bigger

# plot residuals against size
keep=which(allD$Treatment=="No_shrub") #  & allD$logarea.t0 > -1)
plot(allD$logarea.t0[keep],resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$logarea.t0[keep])) # NEG for POSE, zero for HECO (only 3 plots!), zero for PSSP

# plot residuals against W.POSE
keep=which(allD$Treatment=="No_shrub") 
plot(sqrt(allD$W.POSE[keep]),resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$W.POSE[keep])) 

# see if quad level predictions are influenced by variation in #'s of plants per quad
quadN <- aggregate(rep(1,NROW(allD)),by=list(quad=allD$quad,year=allD$year),
                   FUN=sum)
names(quadN)[3]="Density"
allD <- merge(allD,quadN)
keep=which(allD$Treatment=="No_shrub")
plot(allD$Density[keep],resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$Density[keep]))

# see if quad level predictions are influenced by cover of plants per quad
quadCov <- aggregate(allD$logarea.t0,by=list(quad=allD$quad,year=allD$year),
                   FUN=sum)
names(quadCov)[3]="Cover"
allD <- merge(allD,quadCov)
keep=which(allD$Treatment=="No_shrub")
plot(allD$Cover[keep],resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$Cover[keep]))


### PSSP growth models
# first run PSSPgrowth.r

# play with LMER models

m1.lmer <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
              (logarea.t0|year),data=allD)
summary(m1.lmer)

# add size x treatment
m2.lmer <- lmer(logarea.t1~logarea.t0*Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
              (logarea.t0|year),data=allD)
summary(m2.lmer) # nope

m3.lmer<- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO  + W.POSE + W.PSSP+W.PSSP*Treatment+ W.allcov + W.allpts+
              (1|Group)+(logarea.t0|year),data=allD) 
summary(m3.lmer) # SIGNIFICANT

# plot residuals against size
keep=which(allD$Treatment=="No_shrub") #  & allD$logarea.t0 > -1)
plot(allD$logarea.t0[keep],resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$logarea.t0[keep])) 

# plot residuals against W.PSSP
keep=which(allD$Treatment=="No_shrub") 
plot(sqrt(allD$W.PSSP[keep]),resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$W.PSSP[keep])) 


#### HECO growth

m1.lmer <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
              (1|Group)+(logarea.t0|year),data=allD) 
m2.lmer<- lmer(logarea.t1~logarea.t0*Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
              (1|Group)+(logarea.t0|year),data=allD) 
summary(m1.lmer)
summary(m2.lmer) # nope

m3.lmer<- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO +W.HECO:Treatment + W.POSE + W.PSSP+W.allcov + W.allpts+
              (1|Group)+(logarea.t0|year),data=allD) 
summary(m3.lmer) # nope

# plot residuals against size
keep=which(allD$Treatment=="No_shrub")
plot(allD$logarea.t0[keep],resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$logarea.t0[keep])) # NEG for POSE, zero for HECO (only 3 plots!), zero for PSSP

# plot residuals against W.PSSP
keep=which(allD$Treatment=="No_shrub") 
plot(sqrt(allD$W.PSSP[keep]),resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$W.PSSP[keep])) # NEG for POSE, zero for HECO (only 3 plots!), zero for PSSP


###
### look for trend between grass growth residuals and pretreatment ARTR cover
# first run any growth script

# get ARTR pretreatment cover
covD <- read.csv("C:\\Repos\\ExperimentTests\\removals\\QuadYearCover.csv")
covD <- subset(covD, species=="Artemisia tripartita" & Treatment=="No_shrub" & year==2011)
covD <- covD[,c("quad","cover")]

allD <-merge(allD,covD,all.x=T)

m1.lmer <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
              (logarea.t0|year),data=allD)

keep=which(allD$Treatment=="No_shrub")
plot(allD$cover[keep],resid(m1.lmer)[keep])
summary(lm(resid(m1.lmer)[keep]~allD$cover[keep])) # NEG for POSE, zero for HECO (only 3 plots!), zero for PSSP

###
### Analyze growth residuals in plots*years following big ARTR natural mortality events

# get historical data
root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
infile <- paste0(root,"\\ExperimentTests\\data\\idaho\\speciesData\\ARTR\\ARTR_buf5_dorm2.csv")
tmp=read.csv(infile)
