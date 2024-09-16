
#
# x <- data.frame(replicate = c("rep1", "rep2", "rep3"),
#                 batch = c("batch1", "batch1", "batch2"),
#                 treatment = c("A", "B", "C"),
#                 dose = c(1, 10, 20),
#                 v = 0,0,0)
#
#
# y <- get_prior_pc(x = x, prior_eff_base = c(0, 0.01),
#                   prior_eff_rep = c(-2, 2), prior_eff_batch = c(0, 0.5),
#                   prior_eff_group_mu = c(0, 1),
#                   prior_eff_group_sigma = c(0.1, 0.1),
#                   prior_alpha = 3, model = "M")
#
# y$


rgamma(n = 30, shape = 3, rate = 3/)
