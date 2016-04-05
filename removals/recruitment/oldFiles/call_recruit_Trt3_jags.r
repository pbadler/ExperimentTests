
rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/recruitment",sep="")); # modify as needed 

sppList=c("ARTR","HECO","POSE","PSSP")
outfile="recruit_params_m0.csv"
dataDir1 <- paste(root,"/driversdata/data/idaho/speciesData/",sep="")
dataDir2 <- paste(root,"/driversdata/data/idaho_modern/speciesData/",sep="")
#--------------------------------------------------------

# get old recruitment data
Nspp=length(sppList)
for(i in 1:Nspp){
  infile1=paste(dataDir1,sppList[i],"/recArea.csv",sep="")
  tmpD=read.csv(infile1)
  tmpD=tmpD[,c("quad","year","NRquad","totParea","Group")]
  names(tmpD)[3]=paste("R.",sppList[i],sep="")
  names(tmpD)[4]=paste("cov.",sppList[i],sep="")
  if(i==1){
    D=tmpD
  }else{
    D=merge(D,tmpD,all=T)
  }
}
D[is.na(D)]=0  # replace missing values 
D$year <- D$year+ 1900

# get new recruitment data
Nspp=length(sppList)
for(i in 1:Nspp){
  infile1=paste(dataDir2,sppList[i],"/recArea.csv",sep="")
  tmpD=read.csv(infile1)
  tmpD=tmpD[,c("quad","year","NRquad","totParea","Group")]
  names(tmpD)[3]=paste("R.",sppList[i],sep="")
  names(tmpD)[4]=paste("cov.",sppList[i],sep="")
  if(i==1){
    D2=tmpD
  }else{
    D2=merge(D2,tmpD,all=T)
  }
}
D2[is.na(D2)]=0  # replace missing values 

# combine old and new data
D=rbind(D,D2)
rm(D2)

# merge in treatment data
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D <- merge(D,tmp, all.x=T)

# get rid of precip treatments
D <- subset(D,Treatment!="Irrigation" & Treatment!="Drought")

# clean up removal treatment data
ii <- which(D$year>=2011 & D$Treatment=="No_shrub")
D$cov.ARTR[ii] <- 0
D$R.ARTR[ii] <- NA # don't try to estimate ARTR recruitment in these plots
ii <- which(D$year>=2011 & D$Treatment=="No_grass")
D$cov.HECO[ii] <- 0 ; D$cov.POSE[ii] <- 0 ; D$cov.PSSP[ii] <- 0
D$R.HECO[ii] <- NA ; D$R.POSE[ii] <- NA ; D$R.PSSP[ii] <- NA # don't try to estimate grass recruitment in these plots

# calculate mean cover by group and year
tmpD=subset(D,Treatment=="Control") # only use control plots
tmpD=D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD=aggregate(tmpD[,4:NCOL(tmpD)],by=list("year"=tmpD$year,
  "Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)]=paste("Gcov.",sppList,sep="")

D=merge(D,tmpD,all.x=T)


# set up data objects for bugs  
y=as.matrix(D[,c(paste("R.",sppList,sep=""))])
R.tot=rowSums(y)
parents1=as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100
parents2=as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100
year=as.numeric(as.factor(D$year))
Nyrs=length(unique(D$year))
N=dim(D)[1]
Nspp=length(sppList)
Group=as.numeric(D$Group)
Ngroups=length(unique(Group))
modern.control = ifelse(D$Treatment=="Control" & D$year > 2000, 1, 0)
no.shrub = ifelse(D$Treatment=="No_shrub",1,0)
no.grass = ifelse(D$Treatment=="No_grass",1,0)

# plots
pdf("recruit_data.pdf",height=6,width=8)
par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,3,1))
wts=c(0.6,1,0.65,0.9)
for(i in 1:Nspp){
 plot(parents1[,i],y[,i],xlab="Local parents (% cover)",ylab="Recruits",main=sppList[i],pch=1,col=year)
 trueparents=wts[1]*parents1[,i]+(1-wts[1])*parents2[,i]
 plot(trueparents,y[,i],xlab="Mixed parents (% cover)",ylab="Recruits",main=sppList[i],pch=1,col=year)
}
dev.off()


# fit as negative binomial with random effects in WinBUGS
library(boot)
library(R2WinBUGS)

data=list("N","y","parents1","parents2",
  "year","Nyrs","Nspp","Ngroups","Group","modern.control","no.shrub","no.grass")

inits=list(1)
inits[[1]]=list(intcpt.yr=matrix(1,Nyrs,Nspp),intcpt.mu=rep(1,Nspp),intcpt.tau=rep(1,Nspp),
  intcpt.mod=rep(0,Nspp),intcpt.noshrub=rep(0,Nspp), intcpt.nograss=rep(0,Nspp),
  intcpt.gr=matrix(1,Ngroups,Nspp),g.tau=rep(1,Nspp),
  dd=matrix(0,Nspp,Nspp),theta=rep(1,Nspp)) 
inits[[2]]=list(intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.mu=rep(0,Nspp),intcpt.tau=rep(10,Nspp),
  intcpt.mod=rep(0,Nspp),intcpt.noshrub=rep(0,Nspp), intcpt.nograss=rep(0,Nspp),
  intcpt.gr=matrix(0,Ngroups,Nspp),g.tau=rep(0.1,Nspp),
  dd=matrix(0,Nspp,Nspp),theta=rep(2,Nspp))
  
params=c("intcpt.yr","intcpt.mu","intcpt.tau","intcpt.mod","intcpt.noshrub","intcpt.nograss",
  "intcpt.gr","g.tau","dd","theta","u","lambda") 

# try with jags
library(coda)
library(rjags)

modelFile <- "bugs-Trt3.txt"

# out=bugs(data,inits,params,
#   model.file="bugs-Trt3.txt",
#   n.chains=2,
#   n.iter=20000,
#   n.burnin=10000,
#   #n.iter=40000,
#   #n.burnin=20000,
#   n.thin=10, 
#   debug=T,DIC=T,bugs.directory="c:/WinBUGS14/")  
#   
# tmp=grep("lambda",row.names(out$summary))
# A=row.names(out$summary)[tmp]
# B=out$summary[tmp,1]
# lambda=matrix(NA,dim(y)[1],Nspp)
# C=paste(A,"<-",B,sep="")
# eval(parse(n=length(A),text=C))
# lambda[is.na(lambda)]=0
# par(mfrow=c(2,2))
# for(i in 1:Nspp){
#   plot(y[,i],lambda[,i],xlab="Obs",ylab="Pred",main=sppList[i])
# }
# par(mfrow=c(2,2))
# for(i in 1:Nspp){
#   plot(parents1[,i],lambda[,i],xlab="Obs",ylab="Pred",main=sppList[i])
# }
# 
# write.table(out$summary,outfile,row.names=T,sep=",")
# tmp=paste("DIC",out$DIC,sep=",")
# write.table(tmp,outfile,col.names=F,row.names=F,append=T)
