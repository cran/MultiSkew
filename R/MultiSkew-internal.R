.FisherSkew <-
function(data){
#Computes Fisher' s  measure of skewness,that is 
#the third standardized moment of each variable in the dataset.

#PRELIMINARIES

n<-nrow(data) #number of units
d<-ncol(data) #number of variables
Y<-matrix(nrow=n,ncol=d)#initializes the matrix of standardized variables

x.mean<-colMeans(data) #mean vector
m<-sweep(data,2,x.mean)#centered data

x.sd<-apply(data,2,sd)# standard deviation of the original data matrix
x.sdcorretta<-x.sd*sqrt((n-1)/n)

for(j in 1:d){
Y[,j]<-m[,j]/x.sdcorretta[j]
}

uno<-matrix(c(1),nrow=n,ncol=1)#initialization of the matrix uno
A<-matrix(nrow=2,ncol=d)#converts the final table into a matrix

for(j in 1:d){
 z<-Y[,j]#j-th standardized variable
 
	A[1,j]<-j
 A[2,j]<-round(mean(z^3),digits=4)#Fisher skewness of the i-th variable  
 
	}
			
.AB<-A

return(data.frame(A[1:2,],row.names=c("Variables","Fisher Skewness")))#final table as a dataframe
}
