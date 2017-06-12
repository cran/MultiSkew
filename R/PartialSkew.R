PartialSkew <-
function(data){
#Multivariate skewness, as defined in Mori, Rohatgi e Szekely (1993).
#The output contains three objects: Vector, Scalar and pvalue.
#Vector is the vector-valued skewness introduced by Mori et al(1993).
#Scalar is the squared norm of Vector. 
#pvalue is the probability of observing a value of Scalar greater than the observed one,
#when data are normally distributed.
#The measure has been introduced in the paper
#Mori T.F., Rohatgi V.K. and Szekely G.J. (1993).
#On multivariate skewness and kurtosis. Theory Probab. Appl. 38, 547-551.

#PRELIMINARIES

n<-nrow(data) #number of units
d<-ncol(data) #number of variables
A<-matrix(c(0),nrow=d*d,ncol=d)#initializes the sum of tensor products
Scalar<-c()

concentrationmatrix<-solve(cov(data)*(n-1)/n)#inverse of the covariance matrix
aut<-eigen(concentrationmatrix)#spectrum of concentrationmatrix

VV<-matrix(c(0),nrow=d)# initializes the matrix of the eigenvectors
DD<-diag(sort(aut$values,decreasing=FALSE)) #diagonal matrix of eigenvalues

for(i in ncol(aut$vectors):1){
VV<-cbind(VV,aut$vectors[,i])
}# updates the matrix of eigenvectors

VV<-VV[,2:ncol(VV)] #matrix of eigenvectors

DD.sqrt<-solve(sqrt(DD))#inverse of the diagonal matrix of eigenvalues
AA<-VV%*%DD.sqrt%*%t(VV)
Q<-solve(AA)

x.mean<-colMeans(data) #mean vector
m<-sweep(data,2,x.mean)#centered data

Z<-m%*%Q #standardized data

for(i in 1:n){
A<-A+kronecker(kronecker(Z[i,],t(Z[i,])),Z[i,]) #updates the sum of the tensor products
}

T<-A/n #third standardized moment

##vec operator

I<-diag(x = 1, d, d)#Identity matrix (d x d)
V<-I[,1]

for (i in 2:d){
V<-c(V,(I[,i]))
}

Vector<-t(T)%*%V
Scalar<<-t(Vector)%*%Vector
prova<-pchisq(0.5*n*Scalar/(d+2),df=d)
pvalue<-1-prova
Vector<<-round(Vector,digits=4)
Scalar<<-round(Scalar,digits=4)
pvalue<<-round(pvalue,digits=4)

}
