\name{PM10_2006}
\alias{PM10_2006}
\docType{data}
\title{PM10_2006: dataset
%%   ~~ data name/kind ... ~~
}
\description{The PM10 dataset provides an evaluation of PM10 (particulate matter with an aerodynamic equivalent diameter of up 10 m ) concentrations recorded in Italy during year 2006. The variables, collected from 257 stations, are: average (MEAN) and 50th percentile (MEDIAN) for stations which have valid data with a time coverage of at least 50; 98th percentile (98TH) and maximum value (MAX). Stations are classified by region, province and zone (rural, urban, suburban).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("PM10_2006")}
\format{
  A data frame with 257 observations on the following 5 variables.
  \describe{
    \item{\code{zone}}{a factor with levels \code{R} \code{S} \code{U}}
    \item{\code{mean}}{a numeric vector}
    \item{\code{median}}{a numeric vector}
    \item{\samp{98th}}{a numeric vector}
    \item{\code{max}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{APAT (2007) Environmental data yearbook http:www.apat.gov.it site it-IT APAT Pubblicazioni Annuario_dei_Dati_Ambientali 
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Bartoletti, S. and Loperfido, N. (2010). Modelling Air Pollution Data by the Skew-Normal Distribution. Stochastic Environmental Research & Risk Assessment 24, 513-517.

Christiansen, M. and Loperfido, N. (2014). Improved Approximation of the Sum of Random Vectors by the Skew-Normal Distribution. Journal of Applied Probability 51, 466-482.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(PM10_2006)
## maybe str(PM10_2006) ; plot(PM10_2006) ...
}
\keyword{datasets}
