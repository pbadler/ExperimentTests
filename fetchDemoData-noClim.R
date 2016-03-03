fetchDemoData<-function(doSpp,dataDir){
  
   
    nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 
    survDfile=paste(dataDir,"/speciesdata/",doSpp,"/survD.csv",sep="")
    survD=read.csv(file=survDfile)
    D=survD[survD$allEdge==0,];
    D$year=D$year
    D$logarea=log(D$area)
    D$quad=as.character(D$quad)
    D1=D=D[order(D$X),];
    
    #Drop seedlings
    e <- (D$logarea>0)
    D <- D[e,];
    
    
    #########read data on neighbors
    ringD <- read.csv(paste(dataDir,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
    ringD$year<-ringD$year
    
    # merge D with ringD (D contains fewer rows)
    D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
    D=D[order(D$X),]
    rm(ringD)
    row.names(D) <- NULL  
    
    ## pull out annulus data for the focal species  
    sppCols=which(substr(names(D),1,4)==doSpp); 
    sppData=data.frame(survives=D$survives,age=D$age,ID=D$trackID, year=D$year, logarea=D$logarea, Group=as.factor(D$Group), quad=D$quad, as.matrix(D[,sppCols]));  
    ### Change this so sppData includes the response 
    #colnames(sppData)<-c("survives","age","ID","trackID","year","logarea","Group","quad",colnames(sppData)[(1+nonCompLength.s):length(colnames(sppData))])

  intraD<-sppData  # focal spp and intraspecific neighbors
  intraD<-data.frame(intraD)  
  
  dists <- read.csv(paste(dataDir,"/speciesdata/IdahoDistanceWeights.csv",sep="")); 
  dist_wts<- dists[,paste0(doSpp)]
  
  C <- data.matrix(intraD[,grep(paste(doSpp),names(intraD),value=T)]) #matrix of conspecific areas in the annuli 
  W <- C%*%dist_wts; 
  
   dataS<-data.frame(intraD[,c("survives","year","ID","logarea","age","Group","quad")], W)
    
    dataS$year<-as.numeric(paste(19,dataS$year,sep=""))
    
    dataS$year<-as.factor(dataS$year);
    dataS$Group <- as.factor(dataS$Group);
  return(dataS)
}

