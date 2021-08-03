functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real[] alpha, vector beta, vector gamma, 
    int[] site_type
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int sp1_id[end - start + 1] = y_slice[, 3];
    int sp2_id[end - start + 1] = y_slice[, 4];
    int n[end - start + 1] = y_slice[, 5];
    for (i in 1:(end - start + 1)) {
      lp += binomial_logit_lpmf(
        y[i] | n[i], alpha[site_type[site_id[i]]] + beta[site_id[i]] + 
        gamma[sp1_id[i]] + gamma[sp2_id[i]]
      );
    }
    return lp;
  }
}
data{
  int nb_sites; // Number of sites
  int nb_spp; // Number of species
  int nb_int; // Total number of interactions
  int nb_types; // Number of interaction types
  int nb_groups; // Number of functional groups
  int Y_array[nb_int, 5]; // Response variable, site IDs, species IDs, replicates
  int site_type[nb_sites]; // ID of the interaction type of each sites
  int sp_group[nb_spp]; // ID of the functional group of each species
}
parameters{
  // Model parameters
  real alpha[nb_types];
  vector[nb_sites] zbeta;
  vector[nb_spp] zgamma;
  real<lower=0> sigma_beta[nb_types];
  real<lower=0> sigma_gamma[nb_groups];
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_spp] gamma;
  beta = zbeta .* to_vector(sigma_beta[site_type]);
  for (i in 1:nb_spp)
    gamma[i] = zgamma[i] * sigma_gamma[sp_group[i]];
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  zbeta ~ std_normal();
  zgamma ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, gamma, site_type
  );
}
generated quantities{
  // Compute pointwise link (probability of interaction)
  vector[nb_int] link = inv_logit(
    to_vector(alpha[site_type[Y_array[, 2]]]) + beta[Y_array[, 2]] + gamma[Y_array[, 3]] + 
    gamma[Y_array[, 4]]
  );
  
  // Compute pointwise log-likelihood
  vector[nb_int] log_lik;
  for (i in 1:nb_int) {
    log_lik[i] = binomial_lpmf(Y_array[i, 1] | Y_array[i, 5], link[i]);
  }
}
