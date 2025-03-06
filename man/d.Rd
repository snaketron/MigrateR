\name{d}
\alias{d}
\docType{data}
\title{
Example dataset \code{d}
}
\description{
This dataset \code{d} contains simulated cell migration speed data from an 
imaginery experiment. It includes migration speed values (column \code{v}) 
for different samples (column \code{sample}) treated with chemical compounds 
(column \code{compound}) administered at a dost (column \code{dose}) on
experimental plate with identifier (column \code{plate}).
}
\usage{
data("d", package = "cellmig")
}
\format{
  A data frame with 7200 observations on the following 6 variables.
  \describe{
    \item{\code{v}}{a numeric vector with cell migration speeds}
    \item{\code{sample}}{a character vector}
    \item{\code{compound}}{a character vector}
    \item{\code{dose}}{a character vector}
    \item{\code{group}}{a character vector}
    \item{\code{plate}}{a character vector}
  }
}
\details{
The code used to generate this data is in  inst/script/sim_s.R
}
\source{
The code used to generate this data is in  inst/script/sim_s.R
}
\references{
The code used to generate this data is in  inst/script/sim_s.R
}
\examples{
data(d, package = "cellmig")
str(d)
}
\keyword{datasets}
