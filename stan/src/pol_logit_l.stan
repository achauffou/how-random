functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, vector beta, vector gamma
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int pla_id[end - start + 1] = y_slice[, 3];
    int pol_id[end - start + 1] = y_slice[, 4];
    for (i in 1:(end - start + 1)) {
      lp += bernoulli_logit_lpmf(
        y[i] | beta[site_id[i]] + gamma[pla_id[i]]
      );
    }
    return lp;
  }
}
data{
  int nb_sites; // Number of sites
  int nb_pla; // Number of plants
  int nb_pol; // Number of pollinators
  int nb_int; // Total number of interactions
  int Y_array[nb_int, 4]; // Response variable, site IDs, plant IDs, pollinator IDs
}
parameters{
  // Model parameters
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
  vector[nb_sites] zbeta;
  vector[nb_pla] zgamma;
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_pla] gamma;
  beta = zbeta * sigma_beta;
  gamma = zgamma * sigma_gamma;
}
model{
  // Priors
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  zbeta ~ std_normal();
  zgamma ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, beta, gamma
  );
}
