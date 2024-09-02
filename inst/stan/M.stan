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

  vector [max(r)] var_rep;
  vector [max(b)] var_batch;

  vector [max(g)] eff_group;
  vector [max(s)] eff_z;

  vector [max(g)] var_group;
  vector [max(s)] var_z;

  real <lower=0> eff_sigma;
  real <lower=0> var_sigma;
}

transformed parameters {
  vector<lower=0> [max(s)] mu;
  vector<lower=0> [max(s)] phi;

  vector<lower=0> [max(s)] alpha;
  vector<lower=0> [max(s)] beta;

  vector [max(s)] eff_sample;
  vector [max(s)] var_sample;

  var_sample = var_group[g] + var_sigma * var_z;
  eff_sample = eff_group[g] + eff_sigma * eff_z;

  phi = exp(var_rep[r] + var_batch[b] + var_sample);
  mu =  exp(eff_rep[r] + eff_batch[b] + eff_sample);

  alpha = (mu .* mu) ./ (phi);
  beta = mu ./ phi;
}

model {
  eff_rep ~ normal(-2, 2);
  eff_batch ~ normal(-2, 2);
  eff_group ~ normal(0, 1);
  eff_sigma ~ normal(0, 1);
  eff_z ~ std_normal();

  var_rep ~ normal(-2, 2);
  var_batch ~ normal(-2, 2);
  var_group ~ normal(0, 1);
  var_sigma ~ normal(0, 1);
  var_z ~ std_normal();

  y ~ gamma(alpha[s], beta[s]);
}

generated quantities {
  real y_hat_sample [max(s)];
  real log_lik [N];

  y_hat_sample = gamma_rng(alpha, beta);
  for(i in 1:N) {
    log_lik[i] = gamma_lpdf(y[i] | alpha[s[i]], beta[s[i]]);
  }
}
