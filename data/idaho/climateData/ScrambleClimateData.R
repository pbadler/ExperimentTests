rm(list=ls(all=TRUE)); 
setwd("c:/repos/driversData/data/idaho/ClimateData"); 

scrambleClimateYears = function(X) {
    years=unique(X$year); 
    newYearSeq = sample(years,replace=FALSE); 

    X$newYear = X$year; 
    for(j in 1:length(years)) {
        yearData = which(X$year==years[j]);
        X$newYear[yearData] = newYearSeq[j]
    }    
    dev.new(); plot(X$year,X$newYear); 
    dataSeq = 2*365*X$newYear + 33*X$month + X$day; 
    e = order(dataSeq); X=X[e,]; 

    X$year <- X$newYear;
    X <- subset(X,select=(-newYear)); 
    return(X)
}

graphics.off(); 
X=read.csv("interpClim.csv"); 
fX = scrambleClimateYears(X); 

dev.new(); 
par(mfrow=c(4,1),mar=c(3,3,1,1),mgp=c(2,1,0));
e = c(1:5000)+rpois(1,0.2*nrow(X)); 
plot(X$Tmax[e]); plot(fX$Tmax[e]); 
plot(X$Ppt[e]); plot(fX$Ppt[e]); 

