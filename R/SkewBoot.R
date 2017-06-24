SkewBoot <-
function(data,replicates,units,type){
  #source("MaxSkew.R")
  
  MaxSkew<-NULL
  Scalar<-NULL
  MardiaSkewness<-NULL
 rm("MardiaSkewness")
 
  

#Bootstrap distribution of the chosen measure of skewness.
#The function calls the package MaxSkew 1.1, which needs to be downloaded.
#The number of iterations required by the package MaxSkew is set equal to 5.


##OUTPUT
#histogram of the above mentioned bootstrap distribution
#Pvalue= p-value of the chosen skewness measure
#Vector= vector containing the bootstrap replicates of the chosen skewness measure

##INPUT
#data= data matrix
#replicates= number of bootstrap replicates
#units= number of rows in the data matrices sampled from the original data matrix
#if type="Directional" or "Mardia", units is an integer > the number of variables
#If type="Partial" it is an integer > the number of variables + 1
#type="Directional", "Partial" or "Mardia"

fine<-replicates+1

if(type!="Directional"&&type!="Partial"&&type!="Mardia"){
print("ERROR: type must be either Directional, Partial or Mardia")
}
else{
xxx<-matrix(c(0),ncol=1,nrow=fine)#initializes the vector "matrix"

uno<-matrix(c(1),nrow=units,ncol=1)#initializes the vector "uno"
AA<-matrix(ncol=1, nrow=fine)# initializes the bootstrapped skewnesses
#source('C:/Cinzia/Package/Noemi_da modificare/MultiSkew/MultiSkew/R/FisherSkew.R')
#source('C:/Cinzia/Package/Noemi_da modificare/MultiSkew/MultiSkew/R/PartialSkew.R')
#source('C:/Cinzia/Package/Noemi_da modificare/MultiSkew/MultiSkew/R/SkewMardia.R')



Y<-matrix(nrow=fine,ncol=1)#initializes the vector "Y"

for (b in 1:fine){

xB<-matrix(sample(data,size=ncol(data)*units,replace=TRUE),nrow=units)#sampled matrices
if (type=="Directional"){
proj<-MaxSkew(xB,5,2,FALSE)##applies MaxSkew to xB
xxx[b,1]<-FisherSkew(proj)[2,1]#computes Fisher skewness

}

if(type=="Partial"){#applies PartialSkew to xB
Scalar<-PartialSkew(xB)
xxx[b,1]<-Scalar
}

if(type=="Mardia"){#applies SkewMardia to xB

  SkewMardia(xB)
  Mardiaprova<-MardiaSkewness
  xxx[b,1]<-Mardiaprova

}


xxx.mean<-mean(xxx) #mean vector (value)
m<-xxx-xxx.mean #centered data
xxx.sd<-sd(xxx)#standard deviation of the data xxx
Y[,1]<-m/xxx.sd
z<-Y[,1]#the standardized variable



AA[b,1]<-round(mean(z^3),digits=4)#third standardized cumulant of the sample

}


if(type=="Directional"){
hist(AA[2:fine,1],freq=FALSE,main="Histogram of bootstrapped Directional skewness",xlab="Skewness")
}
if(type=="Partial"){
hist(AA[2:fine,1],freq=FALSE,main="Histogram of bootstrapped Partial skewness",xlab="Skewness")
}
if(type=="Mardia"){
hist(AA[2:fine,1],freq=FALSE,main="Histogram of bootstrapped Mardia skewness",xlab="Skewness")
}

print("Vector")
print(AA[2:fine,1])

##in order to compute the p-value of bootstrap
count<-0
for(i in 2:fine){
if(AA[i,1]>= FisherSkew(data)[2,1])
count<-count+1
}
pvalue.Skew<-(count+1)/fine
print("Pvalue")
print(pvalue.Skew)
}

}
