functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real alpha, vector beta, 
    real lambda_pla, real lambda_pol, vector S_pla, vector S_pol, real mu_pla
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int pla_id[end - start + 1] = y_slice[, 3];
    int pol_id[end - start + 1] = y_slice[, 4];
    int n[end - start + 1] = y_slice[, 5];
    int pla_nat[end - start + 1] = y_slice[, 6];
    for (i in 1:(end - start + 1)) {
      lp += binomial_logit_lpmf(
        y[i] | n[i], alpha + beta[site_id[i]] + 
        lambda_pla * S_pla[start + i - 1] + lambda_pol * S_pol[start + i - 1] + 
        mu_pla * pla_nat[i]
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
  int Y_array[nb_int, 6]; // Response variable, site IDs, plant IDs, pollinator IDs, replicates, origin status
  vector[nb_int] S_pla; // Standardized plants bioclimatic suitabilities
  vector[nb_int] S_pol; // Standardized pollinators bioclimatic suitabilities
}
parameters{
  // Model parameters
  real alpha;
  real lambda_pla;
  real lambda_pol;
  real mu_pla;
  vector[nb_sites] zbeta;
  real<lower=0> sigma_beta;
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  beta = zbeta * sigma_beta;
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  lambda_pla ~ normal(0, 1.3);
  lambda_pla ~ normal(0, 1.3);
  mu_pla ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  zbeta ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, lambda_pla, lambda_pol, S_pla, 
    S_pol, mu_pla
  );
}
generated quantities{
  // Compute pointwise link (probability of interaction)
  vector[nb_int] link = inv_logit(
    alpha + beta[Y_array[, 2]] + lambda_pla * S_pla + lambda_pol * S_pol +
    mu_pla * to_vector(Y_array[, 6])
  );
  
  // Compute pointwise log-likelihood
  vector[nb_int] log_lik;
  for (i in 1:nb_int) {
    log_lik[i] = binomial_lpmf(Y_array[i, 1] | Y_array[i, 5], link[i]);
  }
}
