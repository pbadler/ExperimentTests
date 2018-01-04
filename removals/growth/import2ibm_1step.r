# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  intcpt.trt=rep(0,Nspp),slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),slopeXtrt=rep(0,Nspp),
  nb=matrix(0,Nspp,Nspp+2),sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))

for(i in 1:Nspp){

  infile=paste0("growth/",sppList[i],"_growth.csv")

  Gdata=read.csv(infile)
  # main intercept
  Gpars$intcpt[i]=Gdata$X.Intercept.[1]
  # Treatment effects
  tmp=grep("Treatment",names(Gdata))
  if(max.CI==T) {
    Gpars$intcpt.trt[i] <- Gdata$TreatMaxCI[1] # add max treatment intercept
  }else{
    Gpars$intcpt.trt[i] <-  Gdata[1,tmp[1]] # add treatment intercept
  }
  # Treatment effect on intraspecific density dependence
  tmp=grep("Treatment",names(Gdata))
  if(length(tmp)>1){
    tmp=tmp[2]
    Gpars$slopeXtrt[i]=Gdata[1,tmp[1]]
  }

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
  tmp=grep("W.",(names(Gdata)))
  tmp2=grep("Treatment",names(Gdata))
  if(length(tmp2)>1) tmp=tmp[-which(is.element(tmp,tmp2))] # remove trt X W interaction
  Gpars$nb[i,]=as.numeric(Gdata[1,tmp])
  # variance parameters
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth
grow=function(Gpars,doSpp,doGroup,doYear,sizes,crowding,Treatment){
  if(Treatment!="Control") {
      my.intercept=Gpars$intcpt[doSpp] + Gpars$intcpt.trt[doSpp]
      my.nb=Gpars$nb
      my.nb[doSpp,doSpp] + Gpars$slopeXtrt[doSpp]
  }else{
      my.intercept=Gpars$intcpt[doSpp]
      my.nb=Gpars$nb
      my.nb[doSpp,doSpp] 
  }
  mu=my.intercept+Gpars$intcpt.yr[doYear,doSpp]+Gpars$intcpt.gr[doGroup,doSpp]+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*sizes+ sum(my.nb*crowding) 
  if(is.na(mu)) browser()
  mu
}
