data {
  int<lower=0> N;
  int s [N];
  int g [max(s)];
  int r [max(s)];
  int b [max(s)];
  real prior_eff_base [2];
  real prior_eff_rep [2];
  real prior_eff_batch [2];
  real prior_eff_group_mu [2];
  real prior_eff_group_sigma [2];
  real prior_alpha;
}

parameters {
}

transformed parameters {
}

model {
}

generated quantities {
  vector [max(s)] y;
  real eff_base;
  vector [max(r)] eff_rep;
  vector [max(b)] eff_batch;
  vector [max(g)] eff_group_mu;
  vector [max(s)] eff_z;
  real <lower=0> eff_group_sigma;
  real <lower=0> alpha;
  vector<lower=0> [max(s)] mu;
  vector [max(s)] eff_sample;

  eff_base = normal_rng(prior_eff_base[1], prior_eff_base[2]);
  for(i in 1:max(r)) {
    eff_rep[i] = normal_rng(prior_eff_rep[1], prior_eff_rep[2]);
  }
  for(i in 1:max(b)) {
    eff_batch[i] = normal_rng(prior_eff_batch[1], prior_eff_batch[2]);
  }
  for(i in 1:max(g)) {
    eff_group_mu[i] = normal_rng(prior_eff_group_mu[1], prior_eff_group_mu[2]);
  }
  for(i in 1:max(s)) {
    eff_z[i] = std_normal_rng();
  }
  eff_group_sigma = abs(normal_rng(prior_eff_group_sigma[1], prior_eff_group_sigma[2]));
  if(prior_alpha == 0) {
    alpha = gamma_rng(0.01, 0.01);
  } else {
    alpha = prior_alpha;
  }

  for(i in 1:max(s)) {
    eff_sample[i] = eff_group_mu[g[i]] + eff_group_sigma * eff_z[i];
    mu[i] =  exp(eff_base + eff_rep[r[i]] + eff_batch[b[i]] + eff_sample[i]);
    y[i] = gamma_rng(alpha, alpha ./ mu[i]);
  }
}
