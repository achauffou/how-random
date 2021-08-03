functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real[] alpha, vector beta, 
    vector lambda, vector S1, vector S2, int[] site_type, int[] sp_group
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
        lambda[sp_group[sp1_id[i]]] * S1[start + i - 1] +
        lambda[sp_group[sp2_id[i]]] * S2[start + i - 1]
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
  vector[nb_int] S1; // Standardized bioclimatic suitabilities of first species
  vector[nb_int] S2; // Standardized bioclimatic suitabilities of first species
}
parameters{
  // Model parameters
  real alpha[nb_types];
  vector[nb_groups] lambda;
  vector[nb_sites] zbeta;
  real<lower=0> sigma_beta[nb_types];
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  beta = zbeta .* to_vector(sigma_beta[site_type]);
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  lambda ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  zbeta ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, lambda, S1, S2, 
    site_type, sp_group
  );
}
generated quantities{
  // Compute pointwise link (probability of interaction)
  vector[nb_int] link = inv_logit(
    to_vector(alpha[site_type[Y_array[, 2]]]) + beta[Y_array[, 2]] + 
    lambda[sp_group[Y_array[, 3]]] .* S1 + lambda[sp_group[Y_array[, 4]]] .* S2
  );
  
  // Compute pointwise log-likelihood
  vector[nb_int] log_lik;
  for (i in 1:nb_int) {
    log_lik[i] = binomial_lpmf(Y_array[i, 1] | Y_array[i, 5], link[i]);
  }
}
