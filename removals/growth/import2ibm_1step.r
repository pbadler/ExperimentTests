# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  intcpt.trt=rep(0,Nspp),slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
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
  Gpars$nb[i,]=as.numeric(Gdata[1,tmp])
  # variance parameters
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth
grow=function(Gpars,doSpp,doGroup,doYear,sizes,crowding,Treatment){
  trt.Effect=Gpars$intcpt.trt[doSpp]
  if(Treatment=="Control") trt.Effect=0
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+Gpars$intcpt.gr[doGroup,doSpp]+trt.Effect+
    (Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*sizes+ sum(Gpars$nb[doSpp,]*crowding) 
  if(is.na(mu)) browser()
  mu
}
