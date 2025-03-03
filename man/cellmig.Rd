% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{MigrateR}
\alias{MigrateR}
\title{Model-based quantification of cell migration speed}
\usage{
MigrateR(x, control = NULL)
}
\value{
A list which contains:
\itemize{
\item f = fit as rstan object
\item x = processed input
\item p = plot with posterior predictive check
\item s = summary of model parameters (means, 95\% credible intervals, etc.)
}

}
\description{
The functions takes a data.frame, x, as its main input. Meanwhile,
the input list control can be used to configure the MCMC procedure
performed by R-package rstan.  The output is a list which contains:
1) f = fit as rstan object; 2) x = processed input; 3) p = plot
with posterior predictive check; 4) s = summary of model parameters
(means, medians, 95\% credible intervals, etc.).
}
\details{
The input x must have cell entries as rows and the following columns:
\itemize{
\item sample = character of sample ID (s1, s2, s3, etc.)
\item plate = character id of experimental plate (p1, p2, p3, etc.)
\item compound = character treatment compound (A, B, C, etc.)
\item dose = character dose e.g. mmol (0, 1, 5, 10, low, high, mid, etc.)
\item v = numeric cell speed

The input list named \code{control} confures the Markov Chain Monte Carlo
sampling of R-package rstan (see example).
}

}
\examples{
data(d, package = "MigrateR")
o <- MigrateR(x = d,
             control = list(mcmc_warmup = 300,
                            mcmc_steps = 600,
                            mcmc_chains = 3,
                            mcmc_cores = 1,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.9,
                            max_treedepth = 10))
head(o)
}
