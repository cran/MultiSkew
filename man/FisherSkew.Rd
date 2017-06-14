\name{FisherSkew}
\alias{FisherSkew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fisher's  measure of skewness
%%  ~~function to do ... ~~
}
\description{Computes Fisher's  measure of skewness, that is 
the third standardized moment of each variable in the dataset
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FisherSkew(data)
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
\value{Dataframe containing Fisher's  measure of skewness of each variable of the dataset
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
data(PM10_2006)
PM10_2006_matrix<-data.matrix(PM10_2006)
FisherSkew(PM10_2006_matrix)
}