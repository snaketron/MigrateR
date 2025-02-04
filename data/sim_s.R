# Simulate data with severe replicate effects
n_reps <- 3
n_cells <- 20
eff_b <- c(0, 0.2)
eff_d <- c(0, -1, -2, -0.5, -0.05)
eff_c <- c(0, 0.1, 0.5, 1, -0.5, -1) # ratio of dose

x <- c()
for(c in 1:length(eff_c)) {
  for(d in 1:length(eff_d)) {
    for(b in 1:length(eff_b)) {
      u <- rnorm(n = n_reps, mean = eff_b[b] + eff_d[d]*eff_c[c], sd = 0.15)
      for(s in 1:n_reps) {
        y <- rgamma(n = n_cells, shape = 3.1, rate = 3.1/exp(u[s]))
        x <- rbind(x, data.frame(v = y,
                                 sample = paste0("s-", s),
                                 compound = paste0("c-", c),
                                 dose = paste0("d-", d),
                                 batch = paste0("b-", b)))
      }
    }
  }
}

ggplot(data = x)+
  facet_wrap(facets = ~compound)+
  geom_sina(aes(x = dose, y = v, col = batch))
d <- x
save(d, file = "data/d_clean.RData")





