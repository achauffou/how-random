functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum(
    int[,] y_slice, int start, int end, real alpha, vector beta, 
    vector gamma_pla, vector gamma_pol, vector lambda, vector nu, real[,] S_pla, 
    real[,] S_pol, real[] D_pla, real[] D_pol
  ) {
    real lp = 0.0;
    int y[end - start + 1] = y_slice[, 1];
    int site_id[end - start + 1] = y_slice[, 2];
    int pla_id[end - start + 1] = y_slice[, 3];
    int pol_id[end - start + 1] = y_slice[, 4];
    for (i in 1:(end - start + 1)) {
      lp += bernoulli_logit_lpmf(
        y[i] | alpha + beta[site_id[i]] + gamma_pla[pla_id[i]] + gamma_pol[pol_id[i]] + 
        lambda[site_id[i]] * S_pla[pla_id[i], site_id[i]] * S_pol[pol_id[i], site_id[i]] + 
        nu[site_id[i]] * D_pla[pla_id[i]] * D_pol[pol_id[i]]
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
  real S_pla[nb_pla, nb_sites]; // Suitability variable
  real S_pol[nb_pol, nb_sites]; // Suitability variable
  real D_pla[nb_pla]; // Degree variable
  real D_pol[nb_pol]; // Degree variable
}
parameters{
  // Model parameters
  real alpha;
  vector[nb_sites] zbeta;
  vector[nb_pla] zgamma_pla;
  vector[nb_pol] zgamma_pol;
  vector[nb_sites] zlambda;
  vector[nb_sites] znu;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma_pla;
  real<lower=0> sigma_gamma_pol;
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_nu;
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] beta;
  vector[nb_pla] gamma_pla;
  vector[nb_pol] gamma_pol;
  vector[nb_sites] lambda;
  vector[nb_sites] nu;
  beta = zbeta * sigma_beta;
  gamma_pla = zgamma_pla * sigma_gamma_pla;
  gamma_pol = zgamma_pol * sigma_gamma_pol;
  lambda = zlambda * sigma_lambda;
  nu = znu * sigma_nu;
}
model{
  // Priors
  alpha ~ normal(0, 1.3);
  sigma_beta ~ exponential(1);
  sigma_gamma_pla ~ exponential(1);
  sigma_gamma_pol ~ exponential(1);
  sigma_lambda ~ exponential(1);
  sigma_nu ~ exponential(1);
  zbeta ~ std_normal();
  zgamma_pla ~ std_normal();
  zgamma_pol ~ std_normal();
  zlambda ~ std_normal();
  znu ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum, Y_array, grainsize, alpha, beta, gamma_pla, gamma_pol, 
    lambda, nu, S_pla, S_pol, D_pla, D_pol
  );
}
