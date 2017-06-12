\name{SkewBoot}
\alias{SkewBoot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap inference for multivariate skewness measures
%%  ~~function to do ... ~~
}
\description{Computes the bootstrap distribution, its histogram and the corresponding p-value of the chosen measure of multivariate skewness (Mardia, Partial or Directional), using a given number of bootstrap replicates.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SkewBoot(data, replicates, units, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix
%%     ~~Describe \code{data} here~~
}
  \item{replicates}{number of bootstrap replicates
%%     ~~Describe \code{replicates} here~~
}
  \item{units}{number of rows in the data matrices sampled from the original data matrix

%%     ~~Describe \code{units} here~~
}
  \item{type}{"Directional", "Partial" or "Mardia". If type is set equal to
"Directional" or "Mardia", units is an integer greater than the number of variables. If type set equal to "Partial", units  is an integer greater than the number of variables + 1

%%     ~~Describe \code{type} here~~
}
}
\details{The function calls the package MaxSkew 1.1, which needs to be downloaded.
The number of iterations required by the package MaxSkew is set equal to 5.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{histogram}{plot of the above mentioned bootstrap distribution}
\item{Pvalue }{p-value of the chosen skewness measure }
\item{Vector }{vector containing the bootstrap replicates of the chosen skewness measure  }
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
##library(MaxSkew)
##data(PM10_2006)
##PM10_2006_matrix<-data.matrix(PM10_2006)
##SkewBoot(PM10_2006_matrix[,2:5], 50, 50, "Directional")
}