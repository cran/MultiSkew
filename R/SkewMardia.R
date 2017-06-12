SkewMardia <-
function(data){
#Multivariate skewness as defined in Mardia (1970),
#that is the sum of squared elements in the third standardized cumulant
#of the data matrix.
#The output contains two objects:
#MardiaSkewness and pvalue.
#MardiaSkewness is the squared norm of the third cumulant of the standardized data.
#pvalue is the probability of observing a value of MardiaSkewness greater than 
#the observed one, when data are normally distributed.
#The input(data) is a data matrix.
#The measure has been introduced in the paper Mardia, K.V. (1970),
#Measures of multivariate skewness and kurtosis with applications.Biometrika 57, 519-530.

#PRELIMINARIES
  
MardiaSkewness<-c()
pvalue<-c()

n<-nrow(data) #number of units
d<-ncol(data) #number of variables
A<-matrix(c(0),nrow=d*d,ncol=d)#initializes the sum of tensor products

concentrationmatrix<-solve(cov(data)*(n-1)/n)#inverse of the covariance matrix
aut<-eigen(concentrationmatrix)#spectrum of concentrationmatrix

VV<-matrix(c(0),nrow=d)#initializes the matrix of eigenvectors
DD<-diag(sort(aut$values,decreasing=FALSE))#diagonal matrix of eigenvalues

for(i in ncol(aut$vectors):1){# updates the matrix of eigenvectors
VV<-cbind(VV,aut$vectors[,i])
}
VV<-VV[,2:ncol(VV)]#matrix of eigenvectors
DD.sqrt<-solve(sqrt(DD)) #inverse of the diagonal matrix of the eigenvalues 
AA<-VV%*%DD.sqrt%*%t(VV)#positive square root of the covariance matrix
Q<-solve(AA)#positive square root of the concentration  matrix

x.mean<-colMeans(data) #mean vector
m<-sweep(data,2,x.mean)#centered data

Z<-m%*%Q #standardized data

for(i in 1:n){
A<-A+kronecker(kronecker(Z[i,],t(Z[i,])),Z[i,]) #updates the sum of tensor products
}

T<-A/n #third standardized moment


##now we compute Mardia...

P<-t(T)%*%T #product of the transposed third standardized moment and the third 
		#standardized moment
MardiaSkewness<<-sum(diag(P))#Mardia

chiMardia<-pchisq((MardiaSkewness*n)/6,df=(d*(d+1)*(d+2))/6)
pvalue<<-1-chiMardia
}
