functions{
  // Partial sum enables within-chain parallel computation of log-likelihood
  real partial_sum_lpmf(
    int[] y_slice, int start, int end, vector alpha, vector beta_pol, 
    vector beta_pla, vector lambda, vector nu, real[,] S_pla, real[,] S_pol, 
    real[] D_pla, real[] D_pol, int[] site_id, int[] pol_id, int[] pla_id
  ) {
    real lp;
    for (i in start:end) {
      lp += bernoulli_logit_lpmf(
        y_slice[i - start + 1] |
        alpha[site_id[i]] + beta_pol[pol_id[i]] + beta_pla[pla_id[i]] + 
        lambda[site_id[i]] * S_pla[pla_id[i], site_id[i]] * S_pol[pol_id[i], site_id[i]] + 
        nu[site_id[i]] * D_pla[pla_id[i]] * D_pol[pol_id[i]]
      );
    }
    return lp;
  }
}
data{
  int nb_pol; // Number of pollinators
  int nb_pla; // Number of plants
  int nb_sites; // Number of sites
  int nb_int; // Total number of interactions
  int Y[nb_int]; // Response variable (presence of pairwise interaction)
  int site_id[nb_int]; // Site IDs
  int pol_id[nb_int]; // Pollinator IDs
  int pla_id[nb_int]; // Plant IDs
  real S_pla[nb_pla, nb_sites]; // Suitability variable
  real S_pol[nb_pol, nb_sites]; // Suitability variable
  real D_pla[nb_pla]; // Degree variable
  real D_pol[nb_pol]; // Degree variable
  // Note: There are ways to optimize this. One of them is by defining the data 
  // a bit differently. For example, defining the suitabilities as Sp[N] and 
  // Spl[M] (and similarly for the degree) would reduce the size of the data 
  // objects and the model would run faster. But let's try it like it is now 
  // and deal with the problems as they come.
}
parameters{
  // Model parameters
  vector[nb_sites] zalpha;
  vector[nb_pol] zbeta_pol;
  vector[nb_pla] zbeta_pla;
  vector[nb_sites] znu;
  vector[nb_sites] zlambda;
  real alpha_bar;
  real beta_pol_bar;
  real beta_pla_bar;
  real nu_bar;
  real lambda_bar;
  real<lower=0> sigma_a;
  real<lower=0> sigma_bp;
  real<lower=0> sigma_bpl;
  real<lower=0> sigma_l;
  real<lower=0> sigma_n;
}
transformed parameters{
  // Non-centered parametrization
  vector[nb_sites] lambda;
  vector[nb_sites] alpha;
  vector[nb_pol] beta_pol;
  vector[nb_pla] beta_pla;
  vector[nb_sites] nu;
  alpha = zalpha * sigma_a + alpha_bar;
  lambda = zlambda * sigma_l + lambda_bar;
  beta_pol = zbeta_pol * sigma_bp + beta_pol_bar;
  beta_pla = zbeta_pla * sigma_bpl + beta_pla_bar;
  nu = znu * sigma_n + nu_bar;
}
model{
  // Priors
  sigma_a ~ exponential(1);
  sigma_n ~ exponential(1);
  sigma_l ~ exponential(1);
  sigma_bp ~ exponential(1);
  sigma_bpl ~ exponential(1);
  alpha_bar ~ normal(0, 1.3);
  beta_pol_bar ~ normal(0, 1.3);
  beta_pla_bar ~ normal(0, 1.3);
  lambda_bar ~ std_normal();
  nu_bar ~ std_normal();
  zalpha ~ std_normal();
  zbeta_pla ~ std_normal();
  zbeta_pol ~ std_normal();
  znu ~ std_normal();
  zlambda ~ std_normal();
  
  // Within-chain parallelization grainsize (1 lets Stan choose it)
  int grainsize = 1;
  
  // Compute log-likelihood sum in parallel
  target += reduce_sum(
    partial_sum_lpmf, Y, grainsize, alpha, beta_pol, beta_pla, lambda, nu, 
    S_pla, S_pol, D_pla, D_pol, site_id, pol_id, pla_id
  );
  // Note: log-likelihood is not saved in the output (as generated quantity) to
  // save some disk space. That should not be a problem as long as the model is
  // tested. It can be added as generated quantity or computed manually later.
}
