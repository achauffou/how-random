functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real alpha, vector beta, vector gamma, 
    vector lambda, vector SS
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int sp1_id[end - start + 1] = y_slice[, 3];
    int sp2_id[end - start + 1] = y_slice[, 4];
    for (i in 1:(end - start + 1)) {
      lp += bernoulli_logit_lpmf(
        y[i] | alpha + beta[site_id[i]] + gamma[sp1_id[i]] + gamma[sp2_id[i]] + 
        lambda[site_id[i]] * SS[start + i - 1]
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
  int Y_array[nb_int, 4]; // Response variable, site IDs, species IDs
  int site_type[nb_sites]; // ID of the interaction type of each sites
  int sp_group[nb_spp]; // ID of the functional group of each species
  vector[nb_int] SS; // Standardized product of bioclimatic suitabilities
}
parameters{
  // Model parameters
  real alpha;
  real lambda_bar[nb_types];
  vector[nb_sites] zbeta;
  vector[nb_spp] zgamma;
  vector[nb_sites] zlambda;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma[nb_groups];
  real<lower=0> sigma_lambda[nb_types];
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_spp] gamma;
  vector[nb_sites] lambda;
  beta = zbeta * sigma_beta;
  for (i in 1:nb_spp)
    gamma[i] = zgamma[i] * sigma_gamma[sp_group[i]];
  for (i in 1:nb_sites)
    lambda[i] = zlambda[i] * sigma_lambda[site_type[i]] + lambda_bar[site_type[i]];
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  lambda_bar ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  sigma_lambda ~ exponential(1);
  zbeta ~ std_normal();
  zgamma ~ std_normal();
  zlambda ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, gamma, lambda, SS
  );
}
generated quantities{
  // Compute pointwise link (probability of interaction)
  vector[nb_int] link = inv_logit(
    alpha + beta[Y_array[, 2]] + gamma[Y_array[, 3]] + 
    gamma[Y_array[, 4]] + lambda[Y_array[, 2]] .* SS
  );
  
  // Compute pointwise log-likelihood
  vector[nb_int] log_lik;
  for (i in 1:nb_int) {
    log_lik[i] = bernoulli_lpmf(Y_array[i, 1] | link[i]);
  }
}
