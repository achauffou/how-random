functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real[] alpha, vector beta, vector zgamma, 
    real[] sigma_gamma, vector lambda, vector SS, int[] site_type
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int sp1_id[end - start + 1] = y_slice[, 3];
    int sp2_id[end - start + 1] = y_slice[, 4];
    int n[end - start + 1] = y_slice[, 5];
    for (i in 1:(end - start + 1)) {
      lp += binomial_logit_lpmf(
        y[i] | n[i], alpha[site_type[i]] + beta[site_id[i]] + 
        sigma_gamma[site_type[site_id[i]]] * zgamma[sp1_id[i]] * zgamma[sp2_id[i]] + 
        lambda[site_type[site_id[i]]] * SS[start + i - 1]
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
  vector[nb_int] SS; // Standardized product of bioclimatic suitabilities
}
parameters{
  // Model parameters
  real alpha[nb_types];
  vector[nb_types] lambda;
  vector[nb_sites] zbeta;
  vector[nb_spp] zgamma;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma[nb_groups];
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  beta = zbeta * sigma_beta;
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  lambda ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  zbeta ~ std_normal();
  zgamma ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, zgamma, sigma_gamma, lambda, 
    SS, site_type
  );
}
generated quantities{
  // Compute pointwise link (probability of interaction)
  vector[nb_int] link = inv_logit(
    alpha[site_type[Y_array[, 2]]] + beta[Y_array[, 2]] + 
    to_vector(sigma_gamma[site_type[Y_array[, 2]]]) .* zgamma[Y_array[, 3]] .* 
    zgamma[Y_array[, 4]] + lambda[site_type[Y_array[, 2]]] .* SS
  );
  
  // Compute pointwise log-likelihood
  vector[nb_int] log_lik;
  for (i in 1:nb_int) {
    log_lik[i] = binomial_lpmf(Y_array[i, 1] | Y_array[i, 5], link[i]);
  }
}
