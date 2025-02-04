data {
  int<lower=0> N; // number of cells
  vector[N] y;    // cell velocity for N cells
  int s [N];      // sample ID
  int g [max(s)]; // group ID: treatment x dose
  int b [max(s)]; // replicate ID
}

parameters {
  real <lower=0> alpha;
  vector [max(b)] eff_batch;
  vector [max(g)] eff_group;
  real <lower=0> sigma_group;
  vector [max(s)] z_sample;
}

transformed parameters {
  vector<lower=0> [max(s)] mu;
  vector [max(s)] eff_sample;

  eff_sample = eff_group[g] + sigma_group * z_sample;
  mu =  exp(eff_sample + eff_batch[b]);
}

model {
  alpha ~ gamma(0.01, 0.01);
  eff_batch ~ normal(-2.5, 1.5);
  eff_group ~ normal(0, 1);
  sigma_group ~ normal(0, 0.5);
  z_sample ~ std_normal();

  y ~ gamma(alpha, alpha ./ mu[s]);
}

generated quantities {
  real y_hat_sample [max(s)];
  real log_lik [N];

  y_hat_sample = gamma_rng(alpha, alpha ./ mu);
  for(i in 1:N) {
    log_lik[i] = gamma_lpdf(y[i] | alpha, alpha ./ mu[s[i]]);
  }
}
