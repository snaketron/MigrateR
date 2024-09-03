data {
  int<lower=0> N; // number of cells
  vector[N] y;    // cell velocity for N cells
  int s [N];      // sample ID
  int g [max(s)]; // group: treatment x dose combination
  int r [max(s)]; // replicate: index values 1, 2 or 3
  int b [max(s)]; // batch: index values 1 or 2
}

parameters {
  vector [max(r)] eff_rep;
  vector [max(b)] eff_batch;
  vector [max(g)] eff_group;
  vector [max(s)] eff_z;
  real <lower=0> eff_sigma;
  real <lower=0> alpha;
}

transformed parameters {
  vector<lower=0> [max(s)] mu;
  vector [max(s)] eff_sample;

  eff_sample = eff_group[g] + eff_sigma * eff_z;

  mu =  exp(eff_rep[r] + eff_batch[b] + eff_sample);
}

model {
  eff_rep ~ normal(-2, 2);
  eff_batch ~ normal(-2, 2);
  eff_group ~ normal(0, 1);
  eff_sigma ~ normal(0, 1);
  eff_z ~ std_normal();

  alpha ~ gamma(0.01, 0.01);

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
