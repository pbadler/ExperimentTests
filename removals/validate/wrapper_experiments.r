rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

# do control plots
quadList <- c(Q1,Q2,Q3,Q4,Q5,Q6 )  # no shrub
sppList <-  c("ARTR","HECO","POSE","PSSP")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  source("validate/ibm_validate_removal.r")
}

# do no grass plots
quadList <- c(Q48,Q49,Q51,Q55,Q57,Q58,Q60,Q62)  # no grass
sppList <- c("ARTR")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  source("validate/ibm_validate_removal.r")
}

# do no shrub plots
quadList <- c(Q47,Q50,Q52,Q53,Q54,Q56,Q59,Q61)  # no shrub
sppList <-  c("HECO","POSE","PSSP")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  source("validate/ibm_validate_removal.r")
}
