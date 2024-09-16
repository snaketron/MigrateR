# .libPaths(new = "/mnt/nfs/simo/rpack/")
# require(rstan)
# require(cellvel)
# source("../../r_utils/Stats.R")
# load("/mnt/nfs/simo/package/cellvel/data/s_n10_var_nobatch.RData")
# u <- cellvel::get_prior_pc(x = x, control = list(mcmc_steps=5000, mcmc_warmup = 100), prior_eff_base = c(-3, 2), prior_eff_rep=c(0, .5), prior_eff_batch = c(0, .5), prior_alpha = 3, prior_eff_group_mu = c(0, .5), prior_eff_group_sigma = c(0, 0.5))
# q <- rstan::extract(object = u, par = "y")$y
# apply(X = q, MARGIN = 2, FUN = get_hdi, hdi_level = 0.95)


get_prior_pc <- function(x, control, model,
                         prior_eff_base = c(-2.0, 2.0),
                         prior_eff_rep = c(0.0, 1.0),
                         prior_eff_batch = c(0.0, 1.0),
                         prior_eff_group_mu = c(0.0, 1.0),
                         prior_eff_group_sigma = c(0.0, 1.0),
                         prior_alpha = 0) {
  message("prior predictive check... \n")

  # check inputs
  x <- process_input(x)

  # check control
  control <- process_control(control_in = control)

  # transform data
  q <- x[, c("s", "g", "r", "b")]
  q <- q[duplicated(q)==F, ]
  q <- q[order(q$s, decreasing = F),]

  M <- stanmodels$prior_M

  # fit model
  fit <- sampling(object = M,
                  data = list(N = nrow(q), s = q$s, g = q$g, r = q$r, b = q$b,
                              prior_eff_base = prior_eff_base,
                              prior_eff_rep = prior_eff_rep,
                              prior_eff_batch = prior_eff_batch,
                              prior_eff_group_mu = prior_eff_group_mu,
                              prior_eff_group_sigma = prior_eff_group_sigma,
                              prior_alpha = prior_alpha),
                  chains = 1,
                  cores = 1,
                  iter = control$mcmc_steps,
                  warmup = control$mcmc_warmup,
                  algorithm = "Fixed_param",)

  return(fit)
}
