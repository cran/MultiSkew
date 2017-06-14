MinSkew <-
function(data,dimension){
#Reduces sample skewness by projecting the data onto appropriate linear subspaces.
#dimension = number of required projections.
#It is comprised between 2 and the number of the variables.
#The output contains two objects: Linear (linear function of the variables)
#and Projections (projected data).
#The method is described in the paper 
#Loperfido, N. (2014). Linear Transformations to Symmetry.
#Journal  of Multivariate Analysis 129, 186-192.

#PRELIMINARIES

n<-nrow(data) #number of units
d<-ncol(data) #number of variables
A<-matrix(c(0),nrow=d*d,ncol=d)#initializes the sum of tensor products

Projections<-matrix()
Linear<-matrix()
values<-as.list(seq(1:(d-1)))
concentrationmatrix<-solve(cov(data))#inverse of the covariance matrix
aut<-eigen(concentrationmatrix)#spectrum of concentrationmatrix

VV<-matrix(c(0),nrow=d)# initializes the matrix of the eigenvectors

DD<-diag(sort(aut$values,decreasing=FALSE))#diagonal matrix of the eigenvalues

for(i in ncol(aut$vectors):1){
VV<-cbind(VV,aut$vectors[,i])
}# updates the matrix of eigenvectors

VV<-VV[,2:ncol(VV)]#matrix of eigenvectors
DD.sqrt<-solve(sqrt(DD))#inverse of the diagonal matrix of the eigenvalues
AA<-VV%*%DD.sqrt%*%t(VV)
Q<-solve(AA)

x.mean<-colMeans(data) #mean vector
m<-sweep(data,2,x.mean)#centered data

Z<-m%*%Q #standardized data

for(i in 1:n){
A<-A+kronecker(kronecker(Z[i,],t(Z[i,])),Z[i,]) #updates the sum of the tensor products
}
M<-A/n #third standardized moment

autvettval<-eigen(t(M)%*%M)#spectrum of the M matrix (the third standardized moment)

autovettori.M<-matrix(c(0),nrow=d)
autovalori.D<-diag(sort(autvettval$values,decreasing=FALSE))

for(i in ncol(autvettval$vectors):1){
autovettori.M<-cbind(autovettori.M,autvettval$vectors[,i])
}

A<-autovettori.M[,1:dimension+1]
QQ<-Z%*%A 

QQQ<-cbind(-QQ[,1],QQ[,2:dimension])
QQQ<-round(QQQ,digits=4)

Projections<<-QQQ
#Projections<<-round(Projections,digits=4)
Linear<<-A
#print("Linear")
#print(Linear)
#print("Projections")
#print(Projections)
}
