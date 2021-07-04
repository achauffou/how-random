functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, vector beta, vector gamma_pla, 
    vector gamma_pol
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int pla_id[end - start + 1] = y_slice[, 3];
    int pol_id[end - start + 1] = y_slice[, 4];
    int n[end - start + 1] = y_slice[, 5];
    for (i in 1:(end - start + 1)) {
      lp += binomial_logit_lpmf(
        y[i] | n[i], beta[site_id[i]] + gamma_pla[pla_id[i]] + gamma_pol[pol_id[i]]
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
  int Y_array[nb_int, 5]; // Response variable, site IDs, plant IDs, pollinator IDs, replicates
}
parameters{
  // Model parameters
  vector[nb_sites] zbeta;
  vector[nb_pla] zgamma_pla;
  vector[nb_pol] zgamma_pol;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma_pla;
  real<lower=0> sigma_gamma_pol;
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_pla] gamma_pla;
  vector[nb_pol] gamma_pol;
  beta = zbeta * sigma_beta;
  gamma_pla = zgamma_pla * sigma_gamma_pla;
  gamma_pol = zgamma_pol * sigma_gamma_pol;
}
model{
  // Priors
  sigma_beta ~ exponential(1);
  sigma_gamma_pla ~ exponential(1);
  sigma_gamma_pol ~ exponential(1);
  zbeta ~ std_normal();
  zgamma_pla ~ std_normal();
  zgamma_pol ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, beta, gamma_pla, gamma_pol
  );
}
