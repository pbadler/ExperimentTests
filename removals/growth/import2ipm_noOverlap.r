# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),slopeXtrt=rep(0,Nspp),
  nb=matrix(0,Nspp,Nspp),sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))
for(i in 1:Nspp){

  infile=paste0("growth/",sppList[i],"_growth.csv")

  Gdata=read.csv(infile)
  # main intercept
  Gpars$intcpt[i]=Gdata$X.Intercept.[1]
  # random year effects on intercept
  Gpars$intcpt.yr[,i]=Gdata$Intercept.yr
  # group effects
  tmp=which(names(Gdata)=="Group")
  if(length(tmp)>0) Gpars$intcpt.gr[,i]=Gdata$Group[!is.na(Gdata$Group)] 
  # size slope
  Gpars$slope[i]=Gdata$logarea.t0[1]
  # random effects on slope
  tmp=which(names(Gdata)=="logarea.yr")
  if(length(tmp)>0) Gpars$slope.yr[,i]=Gdata[,tmp]
  # get competition coefficients
  tmp=paste("W",sppList,sep=".")
  tmp=which(is.element(names(Gdata),tmp))
  if(length(tmp)>0) Gpars$nb[i,]=as.numeric(Gdata[1,tmp])
  # Treatment effect on intercept
  if(trtEffects==T){
    tmp=grep("Treatment",names(Gdata))[1]
    if(max.CI==T) {
      Gpars$intcpt[i]=Gpars$intcpt[i] + Gdata$TreatMaxCI[1] # add max treatment intercept
    }else{
      Gpars$intcpt[i]=Gpars$intcpt[i] + Gdata[1,tmp[1]] # add treatment intercept
    }
   #Treatment effect on intraspecific density dependence
   tmp=grep("Treatment",names(Gdata))
   if(length(tmp)>1){
    tmp=tmp[2]
    Gpars$nb[i,i]=Gpars$nb[i,i]+Gdata[1,tmp[1]]
   }
  }
  # variance parameters
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth function
G=function(v,u,W,Gpars,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%Gpars$nb[doSpp,]
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}

