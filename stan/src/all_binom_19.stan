functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real[] alpha, vector beta, vector gamma, 
    vector lambda, vector S1, vector S2, int[] site_type, int[] sp_group, vector mu
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int sp1_id[end - start + 1] = y_slice[, 3];
    int sp2_id[end - start + 1] = y_slice[, 4];
    int n[end - start + 1] = y_slice[, 5];
    int sp1_nat[end - start + 1] = y_slice[, 6];
    for (i in 1:(end - start + 1)) {
      lp += binomial_logit_lpmf(
        y[i] | n[i], alpha[site_type[site_id[i]]] + beta[site_id[i]] + 
        gamma[sp1_id[i]] + gamma[sp2_id[i]] + 
        lambda[sp_group[sp1_id[i]]] * S1[start + i - 1] +
        lambda[sp_group[sp2_id[i]]] * S2[start + i - 1] +
        sp1_nat[i] * mu[sp_group[sp1_id[i]]]
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
  int Y_array[nb_int, 6]; // Response variable, site IDs, species IDs, replicates, origin status
  int site_type[nb_sites]; // ID of the interaction type of each sites
  int sp_group[nb_spp]; // ID of the functional group of each species
  vector[nb_int] S1; // Standardized bioclimatic suitabilities of first species
  vector[nb_int] S2; // Standardized bioclimatic suitabilities of first species
}
parameters{
  // Model parameters
  real alpha[nb_types];
  vector[nb_groups] lambda;
  vector[nb_types] mu;
  vector[nb_sites] zbeta;
  vector[nb_spp] zgamma;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma[nb_groups];
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_spp] gamma;
  beta = zbeta * sigma_beta;
  for (i in 1:nb_spp)
    gamma[i] = zgamma[i] * sigma_gamma[sp_group[i]];
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  lambda ~ normal(0, 1.3);
  mu ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  zbeta ~ std_normal();
  zgamma ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, gamma, lambda, S1, S2, 
    site_type, sp_group, mu
  );
}
generated quantities{
  // Compute pointwise link (probability of interaction)
  vector[nb_int] link = inv_logit(
    to_vector(alpha[site_type[Y_array[, 2]]]) + beta[Y_array[, 2]] + gamma[Y_array[, 3]] + 
    gamma[Y_array[, 4]] + lambda[sp_group[Y_array[, 3]]] .* S1 + 
    lambda[sp_group[Y_array[, 4]]] .* S2 + 
    mu[sp_group[Y_array[, 3]]] .* to_vector(Y_array[, 6])
  );
  
  // Compute pointwise log-likelihood
  vector[nb_int] log_lik;
  for (i in 1:nb_int) {
    log_lik[i] = binomial_lpmf(Y_array[i, 1] | Y_array[i, 5], link[i]);
  }
}
