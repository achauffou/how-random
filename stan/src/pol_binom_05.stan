functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real alpha, vector beta, vector lambda, 
    vector SS
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    for (i in 1:(end - start + 1)) {
      lp += bernoulli_logit_lpmf(
        y[i] | alpha + beta[site_id[i]] + lambda[site_id[i]] * SS[start + i - 1]
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
  vector[nb_int] SS; // Standardized product of bioclimatic suitabilities
}
parameters{
  // Model parameters
  real alpha;
  real lambda_bar;
  vector[nb_sites] zbeta;
  vector[nb_sites] zlambda;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_lambda;
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_sites] lambda;
  beta = zbeta * sigma_beta;
  lambda = zlambda * sigma_lambda + lambda_bar;
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  lambda_bar ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  sigma_lambda ~ exponential(1);
  zbeta ~ std_normal();
  zlambda ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, lambda, SS
  );
}
