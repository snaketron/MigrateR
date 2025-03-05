\name{sim}
\alias{sim}
\title{Simulate Data from a Stan Model}
\description{
  Simulates well-level data based on a hierarchical model using a Stan model.
}
\usage{
  sim(N_well_cells, 
      N_plate, 
      N_group, 
      N_well_reps, 
      shape, 
      sigma_bplate, 
      sigma_wplate, 
      alpha_plate, 
      mu_group)
}
\arguments{
  \item{N_well_cells}{Number of cells in each well.}
  \item{N_plate}{Number of plates in the simulation.}
  \item{N_group}{Number of experimental groups (e.g. compounds 
  administrated at different doses).}
  \item{N_well_reps}{Number of technical replicates (wells) per 
  group per plate.}
  \item{shape}{Shape parameter for the Gamma distribution.}
  \item{sigma_bplate}{Between-plate variability in the group 
  effects on cell migration.}
  \item{sigma_wplate}{Within-plate (well) variability in the 
  group effects on cell migration.}
  \item{alpha_plate}{Vector of plate-level intercept effects 
  (baseline effects) on cell migration.}
  \item{mu_group}{Vector of group-level mean effects on cell 
  migration.}
}
\details{
  This function constructs a hierarchical dataset by simulating values 
  using a predefined Stan model. It generates metadata for each well and 
  simulates data using the `sampling` function from the `rstan` package.
}
\value{
  A data frame containing:
  \itemize{
    \item{iteration}{Simulation iteration index.}
    \item{well_id}{Unique well identifier.}
    \item{y}{Simulated response value.}
    \item{group_id}{Experimental group identifier.}
    \item{plate_id}{Plate identifier.}
  }
}
\examples{
ds <- sim(N_well_cells = 100, 
          N_plate = 3, 
          N_group = 5, 
          N_well_reps = 3, 
          shape = 4, 
          sigma_bplate = 0.3, 
          sigma_wplate = 0.1, 
          alpha_plate = c(1, 1.2, 1), 
          mu_group = c(-2, -1, 0, 1, 2))
}
\seealso{
  \code{sampling} from the \pkg{rstan} package.
}
