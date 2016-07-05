
# call from "removal_analysis_wrapper.r"

#============================================================
# SIMULATION PARAMETERS
#============================================================
A=10000 #Area of 100cm x 100cm quadrat
sppList=c("ARTR","HECO","POSE","PSSP")
bigM=c(75,75,50,50)     #Set matrix dimension for each species
maxSize=c(3000,202,260,225)    # in cm^2: PSSP=225 HECO=202  POSE=260  ARTR=3000  # minSize=0.2  cm^2
Nyrs=30
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group

#============================================================
# LOAD VITAL RATE PARAMETERS & FUNCTIONS
#============================================================
Nspp=length(sppList)
# set up survival parameters and function
source("survival/import2ipm_noOverlap.r")
# set up growth parameters and function
source("growth/import2ipm_noOverlap.r")
# set up recruitment parameters and function
source("recruitment/import2ipm.r")

# model spatial group variation (or not)
if(!is.na(doGroup)){
  Spars$intcpt=Spars$intcpt+Spars$intcpt.gr[doGroup,]
  Gpars$intcpt=Gpars$intcpt+Gpars$intcpt.gr[doGroup,]
  Rpars$intcpt.yr=Rpars$intcpt.yr+matrix(Rpars$intcpt.gr[doGroup,],Nyrs,Nspp,byrow=T)
}

#============================================================================================#
# (II) Simulation length, Matrix size and initial vectors
#============================================================================================#

v=v.r=b.r=expv=Cr=Wmat=list(4)
h=r.L=r.U=Ctot=numeric(4)
for(i in 1:Nspp){
  
  # minimum (0.9*minimum size from data) and maximum sizes (1.5*maximum size from data)
  L=log(0.2)
  U=log(maxSize[i]*2)     
  
  # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
  b = L+c(0:bigM[i])*(U-L)/bigM[i] 
  
  # v calculates the middle of each n-equal-sized portion.
  v[[i]] = 0.5*(b[1:bigM[i]]+b[2:(bigM[i]+1)])
  
  # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
  h[i] = v[[i]][2]-v[[i]][1]  
  
  # variables for Wr approximation        Q
  b.r[[i]]=sqrt(exp(b)/pi)
  v.r[[i]]=sqrt(exp(v[[i]])/pi)
  expv[[i]]=exp(v[[i]])
  r.L[i] = sqrt(exp(L)/pi); 
  r.U[i] = sqrt(exp(U)/pi); 
  Wmat[[i]]=matrix(NA,length(v.r[[i]]),Nspp)  # storage of size-specific W values for each focal species
  
} # next species
tmp=range(v.r)
size.range=seq(tmp[1],tmp[2],length=50) # range across all possible sizes

#============================================================================================#
# (III) Utility functions
#============================================================================================#

# load the necessary libraries
library(boot)
library(mvtnorm)
library(msm)
library(statmod)  

## combined kernel
make.K.values=function(v,u,muW, #state variables
  Rpars,rpa,Gpars,Spars,doYear,doSpp)  #growth arguments
{
  f(v,u,Rpars,rpa,doSpp)+S(u,muW,Spars,doYear,doSpp)*G(v,u,muW,Gpars,doYear,doSpp) 
}

# Function to make iteration matrix based only on mean crowding
make.K.matrix=function(v,muW,Rpars,rpa,Gpars,Spars,doYear,doSpp) {
       muW=expandW(v,v,muW)
       K.matrix=outer(v,v,make.K.values,muW,Rpars,rpa,Gpars,Spars,doYear,doSpp)
       return(h[doSpp]*K.matrix)
}

# Function to format the W matrix for the outer product
expandW=function(v,u,W){
   if(dim(W)[1]!=length(u)) stop("Check size of W")
   Nspp=dim(W)[2]
   W=as.vector(W)
   W=matrix(W,length(W),ncol=length(v))
   W=as.vector(t(W))
   W=matrix(W,nrow=length(u)*length(v),ncol=Nspp)
   return(W)
}

# Function to calculate size-dependent crowding, assuming no overlap
wrij=function(r,i,j) {
   return(2*pi*integrate(function(z) z*exp(-alpha[i,j]*(z^2))*Cr[[j]](z-r), r, r + r.U[j])$value
     + pi*Ctot[j]*exp(-alpha[i,j]*((r+r.U[j])^2))/alpha[i,j]); 	
}
Wrij=Vectorize(wrij,vectorize.args="r")


# Function to sum total cover of each species
sumCover=function(v,nt,h,A){
   out=lapply(1:Nspp,function(i,v,nt,h,A) h[i]*sum(nt[[i]]*exp(v[[i]]))/A,v=v,nt=nt,h=h,A=A)
   return(unlist(out))
} 

# Function to sum total density of each species
sumN=function(nt,h){
   out=lapply(1:Nspp,function(i,nt,h) h[i]*sum(nt[[i]]),nt=nt,h=h)
   return(unlist(out))
}

# Function to calculate size variance of each species
varN=function(v,nt,h,Xbar,N){
   out=lapply(1:Nspp,function(i,v,nt,h,Xbar,N) h[i]*sum((exp(v[[i]]-Xbar[i])^2)*nt[[i]])/N[i],v=v,nt=nt,h=h,Xbar=Xbar,N=N)
   return(unlist(out))
}  
              
# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(100),...) {
	nx=length(x); ny=length(y); 
	x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
	y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}
