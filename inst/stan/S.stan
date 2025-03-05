data {
  int<lower=0> N_plate;         // number of plates
  int<lower=0> N_group;         // number groups
  int<lower=0> N_well_reps;     // number of technical replicates
  real <lower=0> shape;         // shape par. of gamma dist.
  real <lower=0> sigma_bplate;  // between-plate variability
  real <lower=0> sigma_wplate;  // within-plate variability
  vector [N_plate] alpha_plate; // plate-specific intercept
  vector [N_group] mu_group;    // average effect of groups (compounds/doses)
}

transformed data {
  int<lower=0> N_well = N_group * N_plate * N_well_reps;
}

parameters {}

transformed parameters {}

model {}

generated quantities {
  real y_hat_sample [N_well];
  vector [N_group] mu_plate_group [N_plate];
  vector<lower=0> [N_well] rate;
  vector [N_well] mu_well;
  int well_id;
  
  well_id = 1;
  for(plate_index in 1:N_plate) {
    for(group_index in 1:N_group) {
      mu_plate_group[plate_index][group_index] = normal_rng(mu_group[group_index], sigma_bplate);
      for(w in 1:N_well_reps) {
        mu_well[well_id] = normal_rng(mu_plate_group[plate_index][group_index], sigma_wplate);
        well_id = well_id + 1;
      }
    }
    rate =  exp(mu_well + alpha_plate[plate_index]);
  }
 
  y_hat_sample = gamma_rng(shape, shape ./ rate);
}
