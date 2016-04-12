
# function to format and output parameters of growth models

formatGrowthPars <-function(model,outfile){
  
  # fit variance
  x=fitted(model)
  y=resid(model)^2
  plot(x,y)
  outVar=nls(y~a*exp(b*x),start=list(a=1,b=0))
  
  #year effects
  params=as.data.frame(ranef(model)$year)
  names(params)=paste(names(params),".yr",sep="")
  year=as.numeric(row.names(params))
  params=cbind(year,params)
  #group effects
  tmp=as.data.frame(ranef(model)$Group)
  if(dim(tmp)[1]>0){
    names(tmp)="Group"
    tmp$GrpName=row.names(tmp)
    tmp[(nrow(tmp)+1):dim(params)[1],]=NA
    params=cbind(params,tmp)
  }
  #fixed effects
  tmp=matrix(NA,dim(params)[1],length(fixef(model)))
  colnames(tmp)=names(fixef(model))
  tmp[1,]=fixef(model)
  params=cbind(params,tmp)
  #variance 
  params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
  params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]
  
  write.table(params,outfile,row.names=F,sep=",")
    
}