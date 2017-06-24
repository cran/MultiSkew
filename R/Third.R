Third <-
function(data,type){
#Third multivariate moment of the data matrix
#containing all moments of order three which can be obtained from the variables
#Some general information about the third multivariate moment of both theoretical
#and emprical distributions are reviewed in the paper
#Loperfido, N. (2015). Singular Value Decomposition of the Third Multivariate Moment.
#Linear Algebra and its Applications 473, 202-216.

##INPUT
#data=data matrix
#type="raw" is the third raw moment
#type="central" is the third central moment
#type="standardized" is the third standardized moment

n<-nrow(data) #number of units
d<-ncol(data) #number of variables

A<-matrix(c(0),nrow=d*d,ncol=d)#initializes the sum of tensor products
ThirdMoment<-NULL
rm("ThirdMoment")
##ThirdMoment<<-matrix()

if(type!="raw"&&type!="central"&&type!="standardized"){
print("ERROR: type must be either raw,central or standardized")
}
else{


if(type=="raw"){
Y<-data
}

if(type=="central"){
x.mean<-colMeans(data) #mean vector
Y<-sweep(data,2,x.mean)#centered data
}
if(type=="standardized"){
concentrationmatrix<-solve(cov(data)*(n-1)/n)
aut<-eigen(concentrationmatrix)
VV<-matrix(c(0),nrow=d)
DD<-diag(sort(aut$values,decreasing=FALSE))

for(i in ncol(aut$vectors):1){
VV<-cbind(VV,aut$vectors[,i])
}
VV<-VV[,2:ncol(VV)]
DD.sqrt<-solve(sqrt(DD))
AA<-VV%*%DD.sqrt%*%t(VV)
Q<-solve(AA)

x.mean<-colMeans(data) #mean vector
m<-sweep(data,2,x.mean)#centered data

Y<-m%*%Q #standardized data
}
}
if(type=="raw"|type=="central"|type=="standardized"){

for(i in 1:n){
A<-A+kronecker(kronecker(Y[i,],t(Y[i,])),Y[i,]) #updates the sum of tensor products
}
M<-A/n # third moment
ThirdMoment<<-round(M,digits=4)

print("Third moment")
print(ThirdMoment)

}
}
