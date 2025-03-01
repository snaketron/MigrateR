data {
  int<lower=0> N; // number of cells
  vector[N] y;    // cell velocity for N cells
  int s [N];      // sample ID
  int g [max(s)]; // group ID: treatment x dose
  int b [max(s)]; // replicate ID
}

parameters {
  real <lower=0> alpha;
  vector <lower=0> [max(b)] sigma_tech;
  real <lower=0> sigma_bio;
  vector [max(g)] mu_bio;
  vector [max(s)] z_tech;
  vector [max(g)] z_bio  [max(g)];
  
}

transformed parameters {
  vector [max(b)] mu_tech [max(g)];
  vector<lower=0> [max(s)] mu;
  vector [max(s)] mulog;
  
  for(j in 1:max(g)) {
    mu_tech[j] = mu_bio[j] + sigma_bio * z_bio[j];
  }
  for(i in 1:max(s)) {
    mulog[i] = mu_tech[g[i]][b[i]] + sigma_tech[b[i]] * z_tech[i];
  }
  
  mu = exp(mulog);
}

model {
  alpha ~ gamma(0.01, 0.01);
  mu_bio ~ normal(0, 1);
  sigma_bio ~ normal(0, 0.5);
  sigma_tech ~ normal(0, 0.5);
  
  z_tech ~ std_normal();
  for(j in 1:max(g)) {
    z_bio[j] ~ std_normal();
  }
  z_tech ~ std_normal();
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
