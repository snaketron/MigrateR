# Simulate data with severe replicate effects
n_reps <- 3
n_cells <- 20
alpha_plate <- c(1, 1.05, 1, 1)
mu_compound <- c(0, 0.1, -1.3, -2, +2, +1.3) 
mu_dose <- c(0, 0.3, 0.8, 0.4, 0.1)

sigma_bplate <- 0.1
sigma_wplate <- 0.05

shape <- 4

x <- c()
# well counter
wc <- 1
for(c in 1:length(mu_compound)) {
  for(d in 1:length(mu_dose)) {
    
    mu_group <- mu_compound[c]*mu_dose[d]
    
    for(p in 1:length(alpha_plate)) {
      u <- rnorm(n = n_reps, mean = alpha_plate[p] + mu_group, sd = sigma_bplate)
      
      for(s in 1:n_reps) {
        z <- rnorm(n = 1, mean = u, sd = sigma_wplate)
        
        y <- rgamma(n = n_cells, shape = shape, rate = shape/exp(z))
        x <- rbind(x, data.frame(v = y,
                                 sample = paste0("w", wc),
                                 compound = paste0("c", c),
                                 dose = paste0("d", d),
                                 group = paste0("c", c, "|d", d),
                                 plate = paste0("p", p)))
        wc <- wc+1
      }
    }
  }
}

d <- x
save(d, file = "data/d.RData", compress = TRUE)




# apply cellmig
cellmig_out <- cellmig(x = d,
             control = list(mcmc_warmup = 250,
                            mcmc_steps = 750,
                            mcmc_chains = 3,
                            mcmc_cores = 3,
                            mcmc_algorithm = "NUTS",
                            adapt_delta = 0.8,
                            max_treedepth = 10))

save(cellmig_out, file = "data/cellmig_out.RData", compress = TRUE)
