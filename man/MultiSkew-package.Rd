\name{MultiSkew-package}
\alias{MultiSkew-package}
\alias{MultiSkew}
\docType{package}
\title{MultiSkew
%\packageTitle{MultiSkew}
}
\description{Computes the third multivariate cumulant of either the raw, centered or standardized data. Computes the main measures of multivariate skewness, together with their bootstrap distributions. Finally, computes the least skewed linear projections of the data
%\packageDescription{MultiSkew}
}
\details{Package: MultiSkew

Type: Package

Title: Measures, Tests and Removes Multivariate Skewness

Version: 1.1

Date: 2017-06-13

Author: Cinzia Franceschini, Nicola Loperfido

Maintainer: Cinzia Franceschini <cinziafranceschini@msn.com>



License: GPL-2



%The DESCRIPTION file:
%\packageDESCRIPTION{MultiSkew}
%\packageIndices{MultiSkew}
%~~ An overview of how to use the package, including the most important functions ~~
%}
%\author{Cinzia Franceschini and Nicola Loperfido
%\packageAuthor{MultiSkew}

%Maintainer:  Cinzia Franceschini <cinziafranceschini@msn.com>
%\packageMaintainer{MultiSkew}
}
\references{Bartoletti, S. and Loperfido, N. (2010). Modelling Air Pollution Data by the Skew-Normal Distribution. Stochastic Environmental Research & Risk Assessment 24, 513-517.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics & Probability Letters 83, 93-99.

Loperfido, N. (2014). Linear Transformations to Symmetry. Journal  of Multivariate Analysis 129, 186-192.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179.

Mardia, K.V. (1970). Measures of multivariate skewness and kurtosis with applications. Biometrika 57, 519-530.

Mori T.F., Rohatgi V.K. and Szekely G.J. (1993). On multivariate skewness and kurtosis. Theory Probab. Appl. 38, 547-551.

}

\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
##data(PM10_2006)
##PM10_2006_matrix<-data.matrix(PM10_2006)

##MinSkew(PM10_2006_matrix[,2:5],4)

##PartialSkew(PM10_2006_matrix[,2:5])

##library(MaxSkew)
##SkewBoot(PM10_2006_matrix[,2:5], 50, 50, "Directional")

##SkewMardia(PM10_2006_matrix[,2:5])

##Third(PM10_2006_matrix[,2:5], "raw")
}