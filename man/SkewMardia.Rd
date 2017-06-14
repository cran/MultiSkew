\name{SkewMardia}
\alias{SkewMardia}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate skewness as defined in Mardia (1970)
%%  ~~function to do ... ~~
}
\description{Sum of squared elements in the third standardized cumulant
of the data matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SkewMardia(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix
%%     ~~Describe \code{data} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{MardiaSkewness}{Squared norm of the third cumulant of the standardized data
}
\item{pvalue }{Probability of observing a value of MardiaSkewness greater than 
the observed one, when data are normally distributed.
}
%% ...
}
\references{Mardia, K.V. (1970),
Measures of multivariate skewness and kurtosis with applications.Biometrika 57, 519-530.
%% ~put references to the literature/web site here ~
}
\author{Cinzia Franceschini and Nicola Loperfido
%%  ~~who you are~~
}
\note{The measure has been introduced in Mardia, K.V. (1970)
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(PM10_2006)
PM10_2006_matrix<-data.matrix(PM10_2006)
SkewMardia(PM10_2006_matrix[,2:5])}