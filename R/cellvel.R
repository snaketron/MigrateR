#' Model-based quantification of cell velocity
#'
#' The functions takes a data.frame, x, as its main input. Meanwhile,
#' the input list control can be used to configure the MCMC procedure
#' performed by R-package rstan.  The output is a list which contains:
#' 1) f = fit as rstan object; 2) x = processed input; 3) s = summary
#' of model parameters (means, medians, 95% credible intervals, etc.).
#'
#' The input x must have cell entries as rows and the following columns:
#' * replicate = id of the biological replicate (e.g. rep1, rep2, rep3, ...)
#' * batch = id of the experimental batch (e.g. plate X, plat Y, ...)
#' * compound = character id of the treatment compound
#' * dose = numeric
#' * v = numeric cell speed
#'
#' @return a list
#' @export
#' @examples
#' data(d, package = "cellvel")
#' o <- cellvel(x = d,
#'              control = list(mcmc_warmup = 300,
#'                             mcmc_steps = 600,
#'                             mcmc_chains = 3,
#'                             mcmc_cores = 1,
#'                             mcmc_algorithm = "NUTS",
#'                             adapt_delta = 0.9,
#'                             max_treedepth = 10))
#' head(o)
cellvel <- function(x, control = NULL, model) {

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

  # transform data
  q <- x[, c("s", "g", "r", "b")]
  q <- q[duplicated(q)==F, ]
  q <- q[order(q$s, decreasing = F),]

  if(model == "M") {
    M <- stanmodels$M
  }
  if(model == "Mp") {
    M <- stanmodels$Mp
  }
  if(model == "Ms") {
    M <- stanmodels$Ms
  }
  if(model == "Ms_pro") {
    M <- stanmodels$Ms_pro
  }

  # fit model
  fit <- sampling(object = M,
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
  l <- x[, c("s", "sample", "g", "group", "treatment", "dose",
             "r", "replicate", "b", "batch")]
  l_s <- l[duplicated(l)==FALSE, ]
  l_r <- l[duplicated(l[,c("r","replicate")])==FALSE,
           c("r","replicate")]
  l_p <- l[duplicated(l[,c("b","batch")])==FALSE,
           c("b","batch")]
  l_g <- l[duplicated(l[,c("g","group")])==FALSE,
           c("g","group", "treatment", "dose",
             "r","replicate","b","batch")]

  # par: eff_rep
  eff_rep <- data.frame(summary(f, par = "eff_rep")$summary)
  eff_rep$r <- 1:nrow(eff_rep)
  eff_rep <- merge(x = eff_rep, y = l_r, by = "r", all.x = TRUE)

  # par: eff_batch
  eff_batch <- data.frame(summary(f, par = "eff_batch")$summary)
  eff_batch$b <- 1:nrow(eff_batch)
  eff_batch <- merge(x = eff_batch, y = l_p, by = "b", all.x = TRUE)

  # par: eff_group_mu
  eff_group_mu <- data.frame(summary(f, par = "eff_group_mu")$summary)
  eff_group_mu$g <- 1:nrow(eff_group_mu)
  eff_group_mu <- merge(x = eff_group_mu, y = l_g, by = "g", all.x = TRUE)

  # par: eff_sample
  eff_sample <- data.frame(summary(f, par = "eff_sample")$summary)
  eff_sample$s <- 1:nrow(eff_sample)
  eff_sample <- merge(x = eff_sample, y = l_s, by = "s", all.x = TRUE)

  # par: eff_group_sigma
  eff_group_sigma <- data.frame(summary(f, par = "eff_group_sigma")$summary)

  # par: mu
  mu <- data.frame(summary(f, par = "mu")$summary)
  mu$s <- 1:nrow(mu)
  mu <- merge(x = mu, y = l_s, by = "s", all.x = TRUE)

  # par: y_hat_sample
  yhat <- data.frame(summary(f, par = "y_hat_sample")$summary)
  yhat$s <- 1:nrow(yhat)
  yhat <- merge(x = yhat, y = l_s, by = "s", all.x = TRUE)

  return(list(eff_rep = eff_rep, eff_batch = eff_batch,
              eff_group_mu = eff_group_mu, eff_sample = eff_sample,
              eff_group_sigma = eff_group_sigma, mu = mu, yhat = yhat))
}


get_profiles <- function(x, hc_link = "average", hc_dist = "euclidean",
                         select_ds, select_ts) {
  eg <- x$s$eff_group_mu
  es <- x$s$eff_sample

  if(missing(select_ds)==FALSE) {
    if(any(!select_ds %in% unique(eg$dose))) {
      stop("selected doses not found in data")
    }
    eg <- eg[eg$dose %in% select_ds, ]
    es <- es[es$dose %in% select_ds, ]
  }

  if(missing(select_ts)==FALSE) {
    if(any(!select_ts %in% unique(eg$treatment))) {
      stop("selected doses not found in data")
    }
    eg <- eg[eg$treatment %in% select_ts, ]
    es <- es[es$treatment %in% select_ts, ]
  }

  q <- acast(data = eg, formula = treatment~dose, value.var = "mean")

  # hclust
  hc <- hclust(dist(q, method = hc_dist), method = hc_link)
  ph <- as.phylo(x = hc)

  tree <- ggtree(ph, linetype='solid')+
    geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col = "black")+
    geom_tippoint(size = 2, fill = "white", shape = 21)+
    geom_tiplab(color='black', as_ylab = T, align = TRUE)+
    layout_rectangular()+
    theme_bw(base_size = 10)+
    scale_x_continuous(labels = abs)

  tree <- revts(tree)

  t <- tree$data
  t <- t[order(t$y, decreasing = FALSE), ]
  tips <- t$label[t$isTip==TRUE]

  q <- eg
  q$treatment <- factor(q$treatment, levels = rev(tips))

  g <- ggplot(data = q)+
    facet_grid(treatment~., switch = "y")+
    geom_hline(yintercept = 0, linetype = "dashed", col = "gray")+
    geom_point(aes(x = dose, y = mean))+
    geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.), width = 0)+
    scale_y_continuous(position = "right", breaks = scales::pretty_breaks(n = 5))+
    theme_bw(base_size = 10)+
    theme(strip.text.y = element_text(margin = margin(0.01,0.01,0.01,0.01, "cm")))


  q <- es[es$treatment %in% q$treatment, ]
  q$treatment <- factor(q$treatment, levels = rev(tips))
  g2 <- ggplot(data = q)+
    facet_wrap(facets = treatment~replicate, nrow = length(unique(q$treatment)), switch = "y", scales = "free_y")+
    geom_errorbar(aes(x = dose, y = mean, ymin = X2.5., ymax = X97.5.), width = 0, alpha = 0.5)+
    geom_line(aes(x = dose, y = mean))+
    geom_point(aes(x = dose, y = mean))+
    scale_y_continuous(position = "right", breaks = scales::pretty_breaks(n = 3))+
    theme_bw(base_size = 10)+
    theme(legend.position = "none",
          #axis.text.y = element_blank(),
          strip.text.y = element_text(margin = margin(0.01,0.01,0.01,0.01, "cm")))


  gs <- (tree|g|g2)+
    # plot_layout(widths = c(1, .8, 4))+
    plot_annotation(tag_levels = 'A')

  gs
  return(gs)
}


get_dose_comparison <- function(x) {

  get_dose_pmax <- function(x, s, e) {
    stats <- c()
    d <- s[s$dose == x, ]
    ts <- sort(unique(d$treatment))

    for(i in 1:(length(ts)-1)) {
      for(j in (i+1):length(ts)) {
        y_i <- e[,d$g[d$dose == x & d$treatment == ts[i]]]
        y_j <- e[,d$g[d$dose == x & d$treatment == ts[j]]]
        u <- y_i-y_j
        hdi <- get_hdi(vec = u, hdi_level = 0.95)

        p_stat <- get_pmax(x = u)

        stats <- rbind(stats, data.frame(treatment_i = ts[i],
                                         treatment_j = ts[j],
                                         contrast = paste0(ts[i], '-vs-', ts[j]),
                                         dose = x,
                                         b_L95 = hdi[1],
                                         b_H95 = hdi[2],
                                         pmax = p_stat,
                                         key = paste0(ts[i], '-', ts[j], '-', x)))
      }
    }
    return(stats)
  }

  e <- rstan::extract(object = o$f, par = "eff_group_mu")$eff_group_mu

  b <- lapply(X = unique(o$x$dose), s = o$s$eff_group_mu, e = e, FUN = get_dose_pmax)
  b <- do.call(rbind, b)
  return(b)
}


get_pmax <- function(x) {
  l <- length(x)
  return(2*max(sum(x<0)/l, sum(x>0)/l)-1)
}

# Description:
# Computes HDI for vector vec and hdi_level (e.g. 0.95)
# Taken (and renamed) from "Doing Bayesian Analysis", section 25.2.3 R code
# for computing HDI of a MCMC sample
get_hdi <- function(vec, hdi_level) {
  # Computes highest density interval from a sample of representative values,
  # estimated as shortest credible interval.
  # Arguments:
  # sampleVec
  # is a vector of representative values from a probability distribution.
  # credMass
  # is a scalar between 0 and 1, indicating the mass within the credible
  # interval that is to be estimated.
  # Value:
  # HDIlim is a vector containing the limits of the HDI
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi_level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}
