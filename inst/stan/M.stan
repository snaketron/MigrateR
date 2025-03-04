data {
  int<lower=0> N;               // number of cells
  int<lower=0> N_well;          // number of wells
  int<lower=0> N_plate;         // number of plates
  int<lower=0> N_plate_group;   // number of plate groups of wells
  int<lower=0> N_group;         // number groups
  vector[N] y;                  // cell velocity for N cells
  int well_id [N];              // well ID
  int plate_id [N_well];        // plate ID
  int plate_group_id [N_well];  // plate ID
  int group_id [N_plate_group]; // group ID: treatment x dose
}

parameters {
  real <lower=0> shape;
  vector [N_plate] alpha_plate;
  vector [N_group] mu_group;
  
  real <lower=0> sigma_bplate;
  real <lower=0> sigma_wplate;
  
  vector [N_well] z_1;
  vector [N_plate_group] z_2;
}

transformed parameters {
  vector<lower=0> [N_well] rate;
  vector [N_well] mu_well;
  vector [N_plate_group] mu_plate_group;

  mu_plate_group = mu_group[group_id] + sigma_bplate * z_2;
  mu_well = mu_plate_group[plate_group_id] + sigma_wplate * z_1; 
  rate =  exp(mu_well + alpha_plate[plate_id]);
}

model {
  shape ~ gamma(0.01, 0.01);
  alpha_plate ~ normal(0, 1);
  mu_group ~ normal(0, 1);
  sigma_bplate ~ normal(0, 0.5);
  sigma_wplate ~ normal(0, 0.5);
  z_1 ~ std_normal();
  z_2 ~ std_normal();

  y ~ gamma(shape, shape ./ rate[well_id]);
}

generated quantities {
  real y_hat_sample [N_well];
  real log_lik [N];

  y_hat_sample = gamma_rng(shape, shape ./ rate);
  for(i in 1:N) {
    log_lik[i] = gamma_lpdf(y[i] | shape, shape ./ rate[well_id[i]]);
  }
}
