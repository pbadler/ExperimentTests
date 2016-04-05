
# function to format and output parameters of growth models

formatSurvPars <-function(model,outfile){

  #year effects
  Intercept.yr=model$summary.random$yearID$mean
  year=model$summary.random$year$ID
  logarea.yr=model$summary.random$year$mean
  params=data.frame(year,Intercept.yr,logarea.yr)
  #group effects
  
  #START HERE
  
  tmp=as.data.frame(ranef(model)$Group)
  if(dim(tmp)[1]>0){
    names(tmp)="Group"
    tmp$GrpName=row.names(tmp)
    tmp[(nrow(tmp)+1):dim(params)[1],]=NA
    params=cbind(params,tmp)
  }
  #fixed effects
  tmp=matrix(NA,dim(params)[1],length(fixef(model)))
  colnames(tmp)=model$names.fixed
  tmp[1,]=m0$summary.fixed[,1]
  params=cbind(params,tmp)

  write.table(params,outfile,row.names=F,sep=",")
    
}