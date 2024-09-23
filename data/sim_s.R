# Simulate data with severe replicate effects
eff_r <- c(-4, -2, -0)
eff_b <- c(0)

eff_d <- c(0, -1, -2, -0.5, -0.05)
eff_t <- c(0, 0.1, 0.5, 1, -0.5, -1) # ratio of dose

x <- c()
for(t in 1:length(eff_t)) {
  for(d in 1:length(eff_d)) {
    for(r in 1:length(eff_r)) {
      u <- rnorm(n = 1, mean = eff_d[d]*eff_t[t], sd = 0.01)
      y <- rgamma(n = 20, shape = 3.1, rate = 3.1/exp(eff_r[r]+eff_b[1]+u))
      x <- rbind(x, data.frame(v = y,
                               replicate = paste0("r-", r),
                               treatment = paste0("t-", t),
                               dose = d,
                               batch = "b-1"))
    }
  }
}

ggplot(data = x)+
  facet_wrap(facets = ~treatment)+
  geom_sina(aes(x = as.factor(dose), y = v, col = replicate))
d <- x
save(d, file = "data/d_severe_rep_effects.RData")






# Simulate data with no replicate/batch effects
eff_r <- c(0, 0, 0)
eff_b <- c(0)

eff_d <- c(0, -1, -2, -0.5, -0.05)
eff_t <- c(0, 0.1, 0.5, 1, -0.5, -1) # ratio of dose

x <- c()
for(t in 1:length(eff_t)) {
  for(d in 1:length(eff_d)) {
    for(r in 1:length(eff_r)) {
      u <- rnorm(n = 1, mean = eff_d[d]*eff_t[t], sd = 0.01)
      y <- rgamma(n = 20, shape = 3.1, rate = 3.1/exp(eff_r[r]+eff_b[1]+u))
      x <- rbind(x, data.frame(v = y,
                               replicate = paste0("r-", r),
                               treatment = paste0("t-", t),
                               dose = d,
                               batch = "b-1"))
    }
  }
}

ggplot(data = x)+
  facet_wrap(facets = ~treatment)+
  geom_sina(aes(x = as.factor(dose), y = v, col = replicate))
d <- x
save(d, file = "data/d_clean.RData")
