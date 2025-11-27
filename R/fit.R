#' @title Wrapped Gaussian Process
#' @description
#' Fits wrapped Gaussian Process (WGP) regression model. Appropriate when response of interest
#' is angular.
#' @param x `numeric` matrix of inputs.
#' @param y `numeric` vector of angular outputs. Assumes units are in radians from 0-2pi.
#' @param beta `numeric` vector of beta coefficients. If `NULL` will estimate
#' in sampler.
#' @param theta `numeric` vector of lengthscale hyperparameters. If `NULL` will estimate
#' in sampler.
#' @param tau2 `numeric` value of scale parameter. If `NULL` will estimate
#' in sampler.
#' @param g `numeric` value of nugget parameter. If `NULL` will estimate
#' in sampler. 
#' @param z `numeric` value of unwrapped GP. Must be numeric vector of same length as y. 
#' If `NULL` will estimate in sampler.
#' @param w `numeric` value of wrapping locations. Must be vector no larger than
#' length of y. If `NULL` will estimate in sampler.
#' @param sigma2 `numeric` value of variance. If `NULL` will estimate in sampler.
#' @param nu `numeric` value of degrees of freedom. If `NULL` will estimate in sampler.
#' @param w_start `numeric` value of initial wrapping locations. Must be vector no
#' larger than length of y.
#' @param method `charater` choice of method for WGP inference, 'ess' or 'jona'.
#' See details section for more information on these methods.
#' @param k_min `numeric` initial specification minimum wrapping number. Only used
#' when `method = ess`.
#' @param k_range `numeric` range of wrapping numbers. Only used when `method = jona`.
#' @param nmcmc `numeric` integer specifying number of MCMC iterations of sampler.
#' @param verb `boolean` whether MCMC progress should be output.
#' @param temp `boolean` whether tempering should be used in MCMC. See details 
#' section for more information.
#' @param settings `list` of MCMC settings. Use `mcmc_settings` to specify this argument. 
#' If `NULL` default settings will be used (see function for details on defaults).
#'
#' @export
#'
#' @returns Object of type `wrapgp`.
wrapgp <- function(x, y, beta = NULL, theta = NULL, tau2 = NULL,
                   g = 1e-4, z = NULL, w = NULL, sigma2 = NULL, nu = NULL,
                   w_start = NULL, method = "ess", k_min = -2, k_range = c(-5, 5),
                   nmcmc = 10000, verb = T, temp = F, settings = NULL){
  if(is.null(settings)) settings <- mcmc_settings(x, y)
  return(gibbs_wrapgp(x, y, beta, theta, tau2, g, z, w, sigma2, nu, w_start,
                      method, k_min, k_range, nmcmc, temp, verb, settings))
}

#' @title Hierarchical wrapped Gaussian Process for RFID Localization
#' @description
#' Fits hierarchical wrapped Gaussian process (WGP) model for purposes of RFID localization.
#' @param x `list` of `numeric` matrices of input frequencies corresponding to each test.
#' Tests with missing distances should be placed at the end.
#' @param y `list` of `numeric` vectors of output phase angles corresponding to each test.
#' Tests with missing distances should be placed at the end.
#' Assumes units are in radians from 0-2pi.
#' @param d `numeric` vector of distances for tests.
#' @param betas `numeric` vector of starting values for slopes. Must be same length as
#' x and y.
#' @param g `numeric` value of nugget parameter. If `NULL` will be sampled from its posterior. 
#' @param delta `numeric` starting values of hierarchical GP. Must be of length `length(x)`.
#' @param mu_delta `numeric` value of mean parameter for hierarchical GP. If `NULL` will be estimated.
#' @param theta_delta `numeric` value of lengthscale parameter for hierachical GP. If `NULL` will be estimated.
#' @param tau2_delta `numeric` value of scale parameter for hierarchical GP. If `NULL` will be estimated.
#' @param s2_beta `numeric` value of common slope variance. If `NULL` will be estimated.
#' @param d_miss `numeric` vector of starting values for distances to impute, if applicable.
#' @param delta_miss `numeric` vector of starting values for hierarchical GP corresponding to missing
#' distances, if applicable.
#' @param method `charater` choice of method for WGP inference, 'ess' or 'jona'.
#' See details section for more information on these methods.
#' @param k_min `numeric` initial specification minimum wrapping number. Only used
#' when `method = ess`.
#' @param k_range `numeric` range of wrapping numbers. Only used when `method = jona`.
#' @param nmcmc `numeric` integer specifying number of MCMC iterations of sampler.
#' @param verb `boolean` whether MCMC progress should be output.
#' @param settings `list` of MCMC settings. Use `mcmc_settings` to specify this argument. 
#' If `NULL` default settings will be used (see function for details on defaults).
#'
#' @export
#'
#' @returns List of chains for sampled hierarchical parameters, as well as a list of
#' `wrapgp` fits for individual tests.
wrapgp_hier <- function(x, y, d, betas = NULL, g = 1e-4, delta = NULL, 
                        mu_delta = NULL, theta_delta = NULL, tau2_delta = NULL, 
                        s2_beta = NULL, d_miss = NULL, delta_miss = NULL, 
                        method = "ess", k_min = -2, k_range = c(-5, 5), 
                        nmcmc = 10000, verb = T, settings = NULL){
  if(is.null(settings)) settings <- mcmc_settings(x[[1]], y[[1]])
  return(gibbs_wrapgp_hier(x, y, d, betas, g, delta, mu_delta, 
                           theta_delta, tau2_delta, s2_beta, d_miss, 
                           delta_miss, method, k_min, k_range, nmcmc, verb, settings))
}