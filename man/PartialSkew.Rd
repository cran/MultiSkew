\name{PartialSkew}
\alias{PartialSkew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PartialSkew
%%  ~~function to do ... ~~
}
\description{Multivariate skewness, as defined in Mori, Rohatgi e Szekely (1993).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PartialSkew(data)
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
\item{Vector }{The vector-valued skewness introduced by  Mori et al (1993)}
  \item{Scalar }{The squared norm of Vector}
  \item{pvalue }{The probability of observing a value of Scalar greater than the observed one,
when data are normally distributed}
%% ...
}
\references{Mori T.F., Rohatgi V.K. and Szekely G.J. (1993).
On multivariate skewness and kurtosis. Theory Probab. Appl. 38, 547-551.
%% ~put references to the literature/web site here ~
}
\author{Cinzia Franceschini and Nicola Loperfido
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##data(PM10_2006)
##PM10_2006_matrix<-data.matrix(PM10_2006)
##PartialSkew(PM10_2006_matrix[,2:5])
}