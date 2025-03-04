#' Model-based quantification of cell velocity
#'
#' The functions takes a data.frame, x, as its main input. Meanwhile,
#' the input list control can be used to configure the MCMC procedure
#' performed by R-package rstan.  The output is a list which contains:
#' 1) f = fit as rstan object; 2) x = processed input; 3) s = summary
#' of model parameters (means, medians, 95% credible intervals, etc.).
#'
#' The input x must have cell entries as rows and the following columns:
#' * sample = unique id of the biological replicate (e.g. S1, S2, S3, ...)
#' * batch = id of the experimental batch (e.g. plate X, plat Y, ...)
#' * compound = character id of the treatment compound
#' * dose = numeric
#' * v = numeric cell speed
#'
#' @return a list
#' @export
#' @examples
#' data(d, package = "cellmig")
#' o <- cellmig(x = d,
#'              control = list(mcmc_warmup = 300,
#'                             mcmc_steps = 600,
#'                             mcmc_chains = 3,
#'                             mcmc_cores = 1,
#'                             mcmc_algorithm = "NUTS",
#'                             adapt_delta = 0.9,
#'                             max_treedepth = 10))
#' head(o)
cellmig <- function(x, control = NULL, model) {

  # check inputs
  x <- process_input(x)

  # check control
  control <- process_control(control_in = control)

  # fit model
  f <- get_fit(x = x, control = control, model = model)

  # get summary
  s <- get_summary(x = x, f = f)

  return(list(f = f, x = x, s = s))
}


get_fit <- function(x, control, model) {
  message("model fitting... \n")

  if(model == "M") {
    M <- stanmodels$M
  }

  # fit model
  fit <- sampling(object = M,
                  data = list(y = x$d$sv, 
                              N = x$d$N[1], 
                              N_well = x$d$N_well[1], 
                              N_plate = x$d$N_plate[1],
                              N_plate_group = x$d$N_plate_group[1],
                              N_group = x$d$N_group[1],
                              well_id = x$d$well_id,
                              plate_id = x$map_w$plate_id,
                              plate_group_id = x$map_w$plate_group_id,
                              group_id = x$map_pg$group_id),
                  chains = control$mcmc_chains,
                  cores = control$mcmc_cores,
                  iter = control$mcmc_steps,
                  warmup = control$mcmc_warmup,
                  algorithm = control$mcmc_algorithm,
                  control = list(adapt_delta = control$adapt_delta,
                                 max_treedepth = control$max_treedepth),
                  refresh = 100)

  return(fit)
}


get_summary <- function(x, f) {
  message("computing posterior summaries...\n")

  # get unique meta data
  l <- x[, c("s", "sample", "g", "group", "compound", "dose", "b", "batch")]
  l_s <- l[duplicated(l)==FALSE, ]
  l_b <- l[duplicated(l[,c("b","batch")])==FALSE, c("b", "batch")]
  l_g <- l[duplicated(l[,c("g","group")])==FALSE,
           c("g", "group", "compound", "dose", "b", "batch")]

  # par: eff_batch
  eff_batch <- data.frame(summary(f, par = "eff_batch")$summary)
  eff_batch$b <- 1:nrow(eff_batch)
  eff_batch <- merge(x = eff_batch, y = l_b, by = "b", all.x = TRUE)

  # par: eff_group
  eff_group <- data.frame(summary(f, par = "eff_group")$summary)
  eff_group$g <- 1:nrow(eff_group)
  eff_group <- merge(x = eff_group, y = l_g, by = "g", all.x = TRUE)

  # par: eff_sample
  eff_sample <- data.frame(summary(f, par = "eff_sample")$summary)
  eff_sample$s <- 1:nrow(eff_sample)
  eff_sample <- merge(x = eff_sample, y = l_s, by = "s", all.x = TRUE)

  # par: sigma_group
  sigma_group <- data.frame(summary(f, par = "sigma_group")$summary)

  # par: mu
  mu <- data.frame(summary(f, par = "mu")$summary)
  mu$s <- 1:nrow(mu)
  mu <- merge(x = mu, y = l_s, by = "s", all.x = TRUE)

  # par: y_hat_sample
  yhat <- data.frame(summary(f, par = "y_hat_sample")$summary)
  yhat$s <- 1:nrow(yhat)
  yhat <- merge(x = yhat, y = l_s, by = "s", all.x = TRUE)

  return(list(eff_batch = eff_batch,
              eff_group = eff_group, 
              eff_sample = eff_sample,
              sigma_group = sigma_group, 
              mu = mu, yhat = yhat))
}

