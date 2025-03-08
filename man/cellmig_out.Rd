\name{cellmig_out}
\alias{cellmig_out}
\docType{data}
\title{
Example \code{cellmig} output object \code{cellmig_out}
}
\description{
This is an example output of \code{cellmig} generated from the example
dataset \code{d}.
}
\usage{
data("cellmig_out", package = "cellmig")
}
\format{
  A list with 3 elements.
  \describe{
    \item{\code{f}}{Stan fit object}
    \item{\code{x}}{processed input}
    \item{\code{s}}{parameter summary list}
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
data(cellmig_out, package = "cellmig")
str(cellmig_out)
}
\keyword{datasets}
