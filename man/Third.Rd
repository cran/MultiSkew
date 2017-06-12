\name{Third}
\alias{Third}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Third multivariate moment of a data matrix
%%  ~~function to do ... ~~
}
\description{It contains all moments of order three which can be obtained from the variables.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Third(data, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data matrix
%%     ~~Describe \code{data} here~~
}
  \item{type}{type="raw" is the third raw moment
  
  type="central" is the third central moment
  
  type="standardized" is the third standardized moment
%%     ~~Describe \code{type} here~~
}
}
\details{Some general information about the third multivariate moment of both theoretical
and emprical distributions are reviewed in Loperfido, N. (2015).

%%  ~~ If necessary, more details than the description above ~~
}
\value{Third moment:
all moments of order three which can be obtained from the variables in "data".

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Loperfido, N. (2015). Singular Value Decomposition of the Third Multivariate Moment.
Linear Algebra and its Applications 473, 202-216.

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
Third(PM10_2006_matrix[,2:5], "raw")
}