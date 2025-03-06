\name{cellmig}
\alias{cellmig}
\title{Model-based quantification of cell migration speed}
\usage{
cellmig(x, control = NULL)
}
\arguments{
\item{x}{main input, data.frame with columns:
\itemize{
\item sample = character of sample ID or well ID (s1, s2, s3, etc.)
\item plate = character id of experimental plate (p1, p2, p3, etc.)
\item compound = character treatment compound (A, B, C, etc.)
\item dose = character dose e.g. mmol (0, 1, 5, 10, low, high, mid, etc.)
\item v = numeric cell migration speed
}
}
\item{control}{input to configure the Markov Chain Monte Carlo (MCMC) 
sampling algorithm, list, with elements and default values:
\itemize{
\item mcmc_warmup = 500
\item mcmc_steps = 1500
\item mcmc_chains = 4
\item mcmc_cores = 1
\item mcmc_algorithm = "NUTS",
\item adapt_delta = 0.8
\item max_treedepth = 10
}
}



}
\value{
A list which contains:
\itemize{
\item f = fit as rstan object
\item x = processed input
\item s = summary of model parameters (means, 95\% credible intervals, etc.)
}

}
\description{
The functions takes a data.frame, x, as its main input. Meanwhile,
the input list control can be used to configure the MCMC procedure
performed by R-package rstan. The output is a list which contains:
1) f = fit as rstan object; 2) x = processed input; 3) s = summary 
of model parameters (means, medians, 95\% credible intervals, etc.).
}
\examples{
data(d, package = "cellmig")
o <- cellmig(x = d,
             control = list(mcmc_warmup = 300,
                            mcmc_steps = 600,
                            mcmc_chains = 3,
                            mcmc_cores = 1,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))
head(o)
}
