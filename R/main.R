#' Model-based quantification of cell velocity
#'
#' The functions takes a data.frame, x, as its main input. Meanwhile,
#' the input list control can be used to configure the MCMC procedure
#' performed by R-package rstan.  The output is a list which contains:
#' 1) f = fit as rstan object; 2) x = processed input; 3) p = plot
#' with posterior predictive check; 4) s = summary of model parameters
#' (means, medians, 95% credible intervals, etc.).
#'
#' The input x must have cell entries as rows and the following columns:
#' * replicate = id of the biological replicate (e.g. rep1, rep2, rep3, ...)
#' * batch = id of the experimental batch (e.g. plate X, plat Y, ...)
#' * group = treatment group (e.g. compound x at dose 10mmol, ...)
#' * v = numeric cell speed
#'
#' @return a list
#' @export
#' @examples
#' data(d, package = "cellvel")
#' o <- cellvel(x = d)
#' head(o)
cellvel <- function(x, control = NULL) {

  # check inputs
  x <- process_input(x)

  # check control
  control <- process_control(control_in = control)

  # fit model
  f <- get_fit(x = x, control = control)

  # get summary
  s <- get_summary(x = x, f = f)

  # get ppc plots
  p <- get_ppc(x = x, s = s)

  return(list(f = f, p = p, x = x, s = s))
}


get_fit <- function(x, control) {
  message("model fitting... \n")

  # transform data
  q <- x[, c("s", "g", "r", "b")]
  q <- q[duplicated(q)==F, ]
  q <- q[order(q$s, decreasing = F),]

  # fit model
  fit <- sampling(object = stanmodels$M,
                  data = list(y = x$sv, N = nrow(x), s = x$s,
                              g = q$g, r = q$r, b = q$b),
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
  l <- x[, c("s", "sample", "g", "group", "r", "replicate", "b", "batch")]
  l_s <- l[duplicated(l)==FALSE, ]
  l_r<-l[duplicated(l[,c("r","replicate")])==FALSE,
         c("r","replicate")]
  l_p<-l[duplicated(l[,c("b","batch")])==FALSE,
         c("b","batch")]
  l_g<-l[duplicated(l[,c("g","group","r","replicate","b","batch")])==FALSE,
         c("g","group","r","replicate","b","batch")]


  # par: eff_rep
  eff_rep <- data.frame(summary(f, par = "eff_rep")$summary)
  eff_rep$r <- 1:nrow(eff_rep)
  eff_rep <- merge(x = eff_rep, y = l_r, by = "r", all.x = TRUE)

  # par: var_rep
  # var_rep <- data.frame(summary(f, par = "var_rep")$summary)
  # var_rep$r <- 1:nrow(var_rep)
  # var_rep <- merge(x = var_rep, y = l_r, by = "r", all.x = TRUE)

  # par: eff_batch
  eff_batch <- data.frame(summary(f, par = "eff_batch")$summary)
  eff_batch$b <- 1:nrow(eff_batch)
  eff_batch <- merge(x = eff_batch, y = l_p, by = "b", all.x = TRUE)

  # par: var_batch
  # var_batch <- data.frame(summary(f, par = "var_batch")$summary)
  # var_batch$b <- 1:nrow(var_batch)
  # var_batch <- merge(x = var_batch, y = l_p, by = "b", all.x = TRUE)

  # par: eff_group
  eff_group <- data.frame(summary(f, par = "eff_group")$summary)
  eff_group$g <- 1:nrow(eff_group)
  eff_group <- merge(x = eff_group, y = l_g, by = "g", all.x = TRUE)

  # par: var_group
  # var_group <- data.frame(summary(f, par = "var_group")$summary)
  # var_group$g <- 1:nrow(var_group)
  # var_group <- merge(x = var_group, y = l_g, by = "g", all.x = TRUE)

  # par: eff_sigma, var_sigma
  # sigma <- data.frame(summary(f, par = c("eff_sigma", "var_sigma"))$summary)
  sigma <- data.frame(summary(f, par = c("eff_sigma"))$summary)

  # par: mu
  mu <- data.frame(summary(f, par = "mu")$summary)
  mu$s <- 1:nrow(mu)
  mu <- merge(x = mu, y = l_s, by = "s", all.x = TRUE)

  # par: phi
  # phi <- data.frame(summary(f, par = "phi")$summary)
  # phi$s <- 1:nrow(phi)
  # phi <- merge(x = phi, y = l_s, by = "s", all.x = TRUE)

  # par: y_hat_sample
  yhat <- data.frame(summary(f, par = "y_hat_sample")$summary)
  yhat$s <- 1:nrow(yhat)
  yhat <- merge(x = yhat, y = l_s, by = "s", all.x = TRUE)

  # return(list(eff_rep = eff_rep, eff_batch = eff_batch, eff_group = eff_group,
  #             var_rep = var_rep, var_batch = var_batch, var_group = var_group,
  #             sigma = sigma, mu = mu, phi = phi, yhat = yhat))
  return(list(eff_rep = eff_rep, eff_batch = eff_batch, eff_group = eff_group,
              sigma = sigma, mu = mu, yhat = yhat))
}


get_ppc <- function(x, s) {
  message("computing ppc...\n")

  g <- ggplot()+
    geom_sina(data = x, aes(x = as.factor(s), y = sv))+
    geom_errorbar(data = s$yhat, aes(x = as.factor(s), y = mean,
                                     ymin = X2.5., ymax = X97.5.),
                  col = "orange")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  return(g)
}

