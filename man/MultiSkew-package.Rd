\name{MultiSkew-package}
\alias{MultiSkew-package}
\alias{MultiSkew}
\docType{package}
\title{Measures, Tests and Removes Multivariate Skewness
%\packageTitle{MultiSkew} 
}
\description{Computes the third multivariate cumulant of either the raw, centered or standardized data. Computes the main measures of multivariate skewness, together with their bootstrap distributions. Finally, computes the least skewed linear projections of the data.
%\packageDescription{MultiSkew}
}
\details{Package: MultiSkew
  
Type: Package

Title: Measures, Tests and Removes Multivariate Skewness

Version: 1.0

Date: 2017-05-31

Author: Cinzia Franceschini and Nicola Loperfido

Maintainer: Cinzia Franceschini <cinziafranceschini@msn.com>

Description: The package computes the third multivariate cumulant of either the raw, centered or standardized data. It also computes the main measures of multivariate skewness, together with their bootstrap distributions.Finally, it computes the least skewed linear projections of the data.

License: GPL (>=2)
%\packageDescription{MultiSkew}
\packageIndices{MultiSkew}
}
%\author{
%\packageAuthor{MultiSkew}
%\packageMaintainer{MultiSkew}
%}
\references{Bartoletti, S. and Loperfido, N. (2010). Modelling Air Pollution Data by the Skew-Normal Distribution. Stochastic Environmental Research & Risk Assessment 24, 513-517.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics & Probability Letters 83, 93-99.

Loperfido, N. (2014). Linear Transformations to Symmetry. Journal  of Multivariate Analysis 129, 186-192.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179.

Mardia, K.V. (1970). Measures of multivariate skewness and kurtosis with applications. Biometrika 57, 519-530.

Mori T.F., Rohatgi V.K. and Szekely G.J. (1993). On multivariate skewness and kurtosis. Theory Probab. Appl. 38, 547-551.
}
\keyword{package}

\examples{
##data(PM10_2006)

##PM10_2006_matrix<-data.matrix(PM10_2006) 

##MinSkew(PM10_2006_matrix[,2:5], 4)

##PartialSkew(PM10_2006_matrix[,2:5])

##SkewMardia(PM10_matrix[,2:5])

##Third(PM10_2006_matrix[,2:5], "raw")

##SkewBoot(PM10_2006_matrix[,2:5], 50, 50, "Directional")
}