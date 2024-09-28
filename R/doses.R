
compare_doses <- function(x) {

  get_dose_pmax <- function(x, s, e) {
    stats <- c()
    d <- s[s$dose == x, ]
    ts <- sort(unique(d$treatment))

    for(i in 1:(length(ts)-1)) {
      for(j in (i+1):length(ts)) {
        y_i <- e[,d$g[d$dose == x & d$treatment == ts[i]]]
        y_j <- e[,d$g[d$dose == x & d$treatment == ts[j]]]
        u <- y_i-y_j

        Mu <- mean(u)
        hdi <- get_hdi(vec = u, hdi_level = 0.95)

        p_stat <- get_pmax(x = u)

        stats <- rbind(stats, data.frame(treatment_i = ts[i],
                                         treatment_j = ts[j],
                                         contrast = paste0(ts[i], '-vs-', ts[j]),
                                         dose = x,
                                         M = Mu,
                                         L95 = hdi[1],
                                         H95 = hdi[2],
                                         pmax = p_stat,
                                         key = paste0(ts[i], '-', ts[j], '-', x)))
      }
    }
    return(stats)
  }

  e <- rstan::extract(object = o$f, par = "eff_group_mu")$eff_group_mu

  b <- lapply(X = unique(o$x$dose),
              s = o$s$eff_group_mu,
              e = e,
              FUN = get_dose_pmax)
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
