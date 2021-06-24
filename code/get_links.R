# Miscellaneous function to compute and save links =============================
#' Compute link if necessary and save it to a file
#' 
compute_save_link <- function(data, rstan_fit, link_fun, res_folder, ...) {
  compute_link <- TRUE
  rstan_file <- file.path(res_folder, "rstan-fit.rds")
  link_file <- file.path(res_folder, "link.rds")
  if (file.exists(link_file)) {
    if (file.info(link_file)$ctime > file.info(rstan_file)$ctime) {
      compute_link <- FALSE
      link <- readRDS(link_file)
    }
  }
  if (compute_link == TRUE) {
    link <- link_fun(data, rstan_fit, ...)
    saveRDS(link, link_file)
  }
  link
}


# Link functions for pollination binomial models ===============================
#' Link function of pollination binomial with centered intercepts only
#'
link.pol_binom_01 <- function(data, fit) {
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      beta[x, data$Y_array$site_id] + gamma_pla[x, data$Y_array$pla_id] +
        gamma_pol[x, data$Y_array$pol_id]
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}

#' Link function of pollination binomial with intercepts only
#'
link.pol_binom_02 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        gamma_pla[x, data$Y_array$pla_id] + gamma_pol[x, data$Y_array$pol_id]
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}

#' Link function of pollination binomial with intercepts and slope
#'
link.pol_binom_03 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  lambda <- as.matrix(fit, pars = "lambda")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        gamma_pla[x, data$Y_array$pla_id] + gamma_pol[x, data$Y_array$pol_id] +
        lambda[x, data$Y_array$site_id] * data$SS
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}

#' Link function of pollination binomial with single lambda for all sites
#'
link.pol_binom_04 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  lambda <- as.matrix(fit, pars = "lambda")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        gamma_pla[x, data$Y_array$pla_id] + gamma_pol[x, data$Y_array$pol_id] +
        lambda * data$SS
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}

#' Link function of pollination binomial with alpha, betas and lambdas
#'
link.pol_binom_05 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  lambda <- as.matrix(fit, pars = "lambda")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        lambda[x, data$Y_array$site_id] * data$SS
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}

#' Link function of pollination binomial with alpha, beta and single lambda
#'
link.pol_binom_06 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  lambda <- as.matrix(fit, pars = "lambda")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] + lambda * data$SS
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}
