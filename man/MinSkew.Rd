\name{MinSkew}
\alias{MinSkew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MinSkew
%%  ~~function to do ... ~~
}
\description{Reduces sample skewness by projecting the data onto appropriate linear subspaces
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MinSkew(data, dimension)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix
%%     ~~Describe \code{data} here~~
}
  \item{dimension}{number of required projections
%%     ~~Describe \code{dimension} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Linear }{linear function of the variables}
   \item{Projections }{projected data}
%% ...
}
\references{Loperfido, N. (2014). Linear Transformations to Symmetry.
Journal  of Multivariate Analysis 129, 186-19
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
data(PM10_2006)
PM10_2006_matrix<-data.matrix(PM10_2006)
MinSkew(PM10_2006_matrix[,2:5],4)
}