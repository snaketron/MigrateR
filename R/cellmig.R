cellmig <- function(x, control = NULL) {

  # check inputs
  x <- process_input(x)

  # check control
  control <- process_control(control_in = control)

  # fit model
  f <- get_fit(x = x, control = control)

  # get summary
  s <- get_summary(x = x, f = f)

  return(list(f = f, x = x, s = s))
}

get_fit <- function(x, control) {
  message("model fitting... \n")

  M <- stanmodels$M

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
                  refresh = 500)

  return(fit)
}

get_summary <- function(x, f) {
  message("computing posterior summaries...\n")

  x <- x$d
  
  # get meta data
  l <- x[, c("well_id", "well", "group_id", "group", 
             "compound", "dose", "plate_id", "plate", 
             "plate_group_id", "plate_group")]
  meta_well <- l[duplicated(l)==FALSE, ]
  meta_plate <- l[duplicated(l[, c("plate", "plate_id")])==FALSE, 
                  c("plate", "plate_id")]
  meta_group <- l[duplicated(l[, c("group", "group_id")])==FALSE, 
                  c("group", "group_id", "compound", "dose", 
                    "plate_id", "plate")]
  meta_plate_group <- l[duplicated(l[, c("group", "group_id", 
                                         "plate", "plate_id",
                                         "plate_group", 
                                         "plate_group_id")])==FALSE, 
                        c("well_id", "group", "group_id", "plate", "plate_id",
                          "compound", "dose", "plate_group", "plate_group_id")]
  
  # par: alpha_plate
  alpha_plate <- data.frame(summary(f, par = "alpha_plate")$summary)
  alpha_plate$plate_id <- 1:nrow(alpha_plate)
  alpha_plate <- merge(x = alpha_plate, y = meta_plate, 
                       by = "plate_id", all.x = TRUE)
  
  # par: mu_group
  mu_group <- data.frame(summary(f, par = "mu_group")$summary)
  mu_group$group_id <- 1:nrow(mu_group)
  mu_group <- merge(x = mu_group, y = meta_group, by = "group_id", all.x = TRUE)

  # par: mu_plate_group
  mu_plate_group <- data.frame(summary(f, par = "mu_plate_group")$summary)
  mu_plate_group$plate_group_id <- 1:nrow(mu_plate_group)
  mu_plate_group <- merge(x = mu_plate_group, y = meta_plate_group, 
                          by = "plate_group_id", all.x = TRUE)
  
  # par: mu_well
  mu_well <- data.frame(summary(f, par = "mu_well")$summary)
  mu_well$well_id <- 1:nrow(mu_well)
  mu_well <- merge(x = mu_well, y = meta_well, by = "well_id", all.x = TRUE)

  # par: sigma_bplate
  sigma_bplate <- data.frame(summary(f, par = "sigma_bplate")$summary)
  sigma_wplate <- data.frame(summary(f, par = "sigma_wplate")$summary)
  
  # par: rate
  rate <- data.frame(summary(f, par = "rate")$summary)
  rate$well_id <- 1:nrow(rate)
  rate <- merge(x = rate, y = meta_well, by = "well_id", all.x = TRUE)

  # par: shape
  shape <- data.frame(summary(f, par = "shape")$summary)
  
  # par: y_hat_sample
  yhat <- data.frame(summary(f, par = "y_hat_sample")$summary)
  yhat$well_id <- 1:nrow(yhat)
  yhat <- merge(x = yhat, y = meta_well, by = "well_id", all.x = TRUE)

  return(list(alpha_plate = alpha_plate,
              mu_group = mu_group, 
              mu_plate_group = mu_plate_group,
              mu_well = mu_well,
              sigma_bplate = sigma_bplate,
              sigma_wplate = sigma_wplate,
              rate = rate, 
              shape = shape,
              yhat = yhat))
}
