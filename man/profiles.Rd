\name{get_profiles}
\alias{get_profiles}
\title{Visualization of compound effect profiles (cellmig results)}
\usage{
get_profiles(x,
             hc_link = "average",
             hc_dist = "euclidean",
             select_ds,
             select_ts)
}
\value{
A patchwork plot containing:
\itemize{
  \item \strong{A}: Hierarchical dendrogram comparing compound effect curves.
  \item \strong{B}: Mean effects of compounds and doses on cell migration 
  with 95\% credible intervals.
  \item \strong{C}: Mean effects of compounds and doses on cell migration 
  in each replicate with 95\% credible intervals.
}
}
\description{
The function takes as its main input (\code{x}) the output of the 
\code{cellmig} function. Users can select a subset of compounds using 
\code{select_ts} and doses using \code{select_ds}. To construct the
hierarchical dendrogram, one has to select a distance metric, 
\code{hc_dist}, (e.g. \code{hc_dist} = "euclidean"); and a linkage 
function, \code{hc_link}, (e.g. \code{hc_link} = "average").
}
\details{
The function \code{get_profiles} visualizes compound effect profiles, 
identifying compounds with similar effects on cell migration across 
different doses. These compounds are clustered together in the 
hierarchical dendrogram.
}
\examples{
data(d, package = "cellmig")
#o <- cellmig(x = d,
#             control = list(mcmc_warmup = 300,
#                            mcmc_steps = 1500,
#                            mcmc_chains = 3,
#                            mcmc_cores = 3,
#                            mcmc_algorithm = "NUTS",
#                            adapt_delta = 0.8,
#                            max_treedepth = 10))
#
data(cellmig_out, package = "cellmig")
p <- get_profiles(x = cellmig_out,
                  hc_link = "average",
                  hc_dist = "euclidean")
p
}
