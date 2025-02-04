
compare_doses <- function(x, select_ds, select_cs) {

  e <- rstan::extract(object = x$f, par = "eff_group")$eff_group
  ds <- unique(x$s$eff_group$dose)
  cs <- unique(x$s$eff_group$compound)

  if(missing(select_ds)==FALSE) {
    if(any(!select_ds %in% ds)) {
      stop("selected doses not found in data")
    }
    ds <- ds[ds %in% select_ds]
  }

  if(missing(select_cs)==FALSE) {
    if(any(!select_cs %in% cs)) {
      stop("selected compounds not found in data")
    }
    cs <- cs[cs %in% select_cs]
  }

  s <- x$s$eff_group
  s <- s[s$compound %in% cs & s$dose %in% ds,]

  b <- lapply(X = ds, s = s, e = e, FUN = get_dose_pmax)
  b <- do.call(rbind, b)
  return(b)
}


get_dose_pmax <- function(x, s, e) {
  stats <- c()
  d <- s[s$dose == x, ]
  cs <- sort(unique(d$compound))

  for(i in 1:(length(cs)-1)) {
    for(j in (i+1):length(cs)) {
      y_i <- e[,d$g[d$dose == x & d$compound == cs[i]]]
      y_j <- e[,d$g[d$dose == x & d$compound == cs[j]]]
      u <- y_i-y_j

      Mu <- mean(u)
      hdi <- get_hdi(vec = u, hdi_level = 0.95)

      p_stat <- get_pmax(x = u)

      stats <- rbind(stats, data.frame(compound_i = cs[i],
                                       compound_j = cs[j],
                                       contrast = paste0(cs[i], '-vs-', cs[j]),
                                       dose = x,
                                       M = Mu,
                                       L95 = hdi[1],
                                       H95 = hdi[2],
                                       pmax = p_stat,
                                       key = paste0(cs[i], '-', cs[j], '-', x)))
    }
  }
  return(stats)
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
