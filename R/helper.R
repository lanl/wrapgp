# wrapped residual
wrap_resid <- function(phi1, phi2) min(abs(phi1 - phi2), abs(phi1 + 2*pi - phi2), abs(phi1 - (phi2 + 2*pi)))

# wrapped rmse
#' @title Circular Root-Mean Squared Error
#' @description Calculates root mean squared-error (RMSE) for predictions of angular responses.
#' @param y `numeric` vector of angular responses.
#' @param mu `numeric` vector of predictions.
#' @export
#' @returns `numeric` Circular RMSE value.
rmse_circ <- function(y, mu) sqrt(mean((sapply(1:length(y), function(i) wrap_resid(y[i], mu[i]))^2)))

# wrapped crps
#' @title Circular Continuous Ranked Probability Score
#' @description Calculates continuous ranked probability score (CRPS) based on 
#' posterior predictive samples of angular responses.
#' @param y `numeric` vector of angular responses.
#' @param draws `numeric` matrix of posterior predictive samples for y.
#' @export
#' @returns `numeric` Circular CRPS value.
crps_circ <- function(y, draws){
  for(i in 1:nrow(draws)){
    for(j in 1:ncol(draws)){
      if(abs(draws[i,j] - 2*pi - y[j]) < abs(draws[i,j] - y[j])){
        draws[i,j] <- draws[i,j] - 2*pi
      }else{
        if(abs(draws[i,j] + 2*pi - y[j]) < abs(draws[i,j] - y[j])){
          draws[i,j] <- draws[i,j] + 2*pi
        }else{
          next
        }
      }
    }
  }
  return(mean(scoringRules::crps_sample(y, t(draws))))
}

# thin chains
#' @title MCMC Chain Trimming
#' @description
#' Trims Markov chain Monte Carlo (MCMC) burn-in and thins chains for sampled parameters in wrapped GP fit.
#' @param fit object of class `wrapgp`.
#' @param burn `numeric` positive integer specifying how many initial MCMC
#' samples to discard as burn-in.
#' @param thin `numeric` positive integer specifying thinning in rate.
#' @export
#' @returns a ``wrapgp`` object with trimmed MCMC chains.
trim <- function(fit, burn, thin){
  nmcmc <- nrow(fit@k)
  inds  <- seq(burn + 1, nmcmc, thin)
  # extract mcmc results
  fit@w      <- fit@w[inds]
  fit@k      <- fit@k[inds,,drop = F]
  fit@z      <- fit@z[inds,,drop = F]
  fit@beta   <- fit@beta[inds,,drop = F]
  fit@mu     <- fit@mu[inds,,drop = F]
  fit@theta  <- fit@theta[inds]
  fit@tau2   <- fit@tau2[inds]
  fit@g      <- fit@g[inds]
  fit@ll_y   <- fit@ll_y[inds]
  fit@ll_z   <- fit@ll_z[inds]
  fit@lambda <- fit@lambda[inds]
  fit@sigma2 <- fit@sigma2[inds]
  fit@nu     <- fit@nu[inds]
  return(fit)
}

# estimate slope
#' @title Slope Estimation
#' @description
#' Estimates slope for empirical prior choice for beta inference in WGP.
#' @param x `numeric` matrix of real-valued input.
#' @param y `numeric` vector of angular output.
#' @param num_points number of points to construct each line segment for slope
#' estimation. Default is 10.
#' @export
#' @returns `numeric` slope estimate.
slope_est <- function(x, y, num_points = 10){
  n <- length(x)
  num_lines <- n %/% num_points
  slopes <- c()
  for(i in 1:num_lines){
    inds <- (num_points*(i - 1) + 1):(num_points*i)
    slope <- sum((x[inds] - mean(x[inds])) * (y[inds] - mean(y[inds]))) / 
      sum((x[inds] - mean(x[inds]))^2)
    if(is.nan(slope)) next
    if(slope > 0) slopes <- c(slopes, slope)
  }
  if(length(slopes) == 0){
    return(0)
  }else{
    return(mean(slopes))
  }
}

# mcmc settings
#' @title MCMC Settings
#' @description Creates list of settings for priors and other MCMC settings
#' to pass to `wrapgp` and `wrapgp_hier`.
#' @param x `numeric` matrix of inputs.
#' @param y `numeric` vector of outputs.
#' @param theta list of prior choice for lengthscale ("uniform" or "gamma"), prior parameters a and b, and step-size u for Metropolis-Hastings proposal. 
#' @param g list of prior parameters a and b for nugget, and step-size u for Metropolis-Hastings proposal.
#' @param tau2 list of prior parameters a and b for scale.
#' @param beta list of mean structure ("constant" or "linear"), prior mean mu0, and
#' covariance Sigma0 for beta.
#' @param sigma2 list of prior parameters a and b for error variance, and step-size u for Metropolis-Hastings proposal.
#' @param nu list of prior parameter a for degrees of freedom, and step-size u for Metropolis-Hastings proposal.
#' @param ref `boolean` of whether reference prior is used for scale.
#' @param nug_onepar `boolean` of whether nugget should be scaled by scale or not.
#' @param alpha `numeric` value of correction factor for slope. Should be greater than 0
#' and no more than 1.
#' @param joint_update `boolean` of whether scale and lengthscale should be sample jointly
#' or individually.
#' @param lambda `numeric` value of tempering parameter. Should be greater than 0
#' and no greater than 1.
#' @param m `numeric` number of ladder "rungs" in tempering procedure. Should be positive integer greater
#' than 1.
#' @param lambda_min `numeric` lower bound for lambda. Should be greater than 0
#' and no greater than 1.
#' @param c0 `numeric` hyperparameter for tempering. Should be greater than 0.
#' @param n0 `numeric` hyperparameter for tempering. Should be greater than 0.
#' @param dist `string` specification of Gaussian ("normal") or Student's-t ("student") likelihood assumption.
#' @export
#' @returns list of MCMC settings.
mcmc_settings <- function(x, y, theta = list(), g = list(), tau2 = list(), beta = list(),
                          sigma2 = list(), nu = list(), ref = F, nug_onepar = F, alpha = 0, joint_update = F,
                          lambda = 1, m = 5, lambda_min = 0.1, c0 = 1, n0 = 1,
                          dist = "student"){
  if(is.null(theta$prior)) theta$prior <- "gamma"
  if(is.null(theta$a)) theta$a <- 2.5
  if(is.null(theta$b)) theta$b <- 1.5
  if(is.null(theta$u)) theta$u <- 2

  if(is.null(g$prior)) g$prior <- "gamma"
  if(is.null(g$a)) g$a <- 1.1
  if(is.null(g$b)) g$b <- 10
  if(is.null(g$u)) g$u <- 2

  if(is.null(tau2$prior)) tau2$prior <- "igamma"
  if(is.null(tau2$a)) tau2$a <- 1
  if(is.null(tau2$b)) tau2$b <- 1

  if(is.null(beta$mean_struct)) beta$mean_struct <- "linear"
  if(is.null(beta$mu0)){
    if(beta$mean_struct == "constant"){
      beta$mu0 <- 0
    }else{
      beta$mu0 <- c(0, slope_est(x, y))
    }
  }
  if(is.null(beta$Sigma0)){
    if(beta$mean_struct == "constant"){
      beta$Sigma0 <- 1
    }else{
      beta$Sigma0 <- diag(c(1, 10))
    }
  }
  if(is.null(beta$mu_s2a)) beta$mu_s2a <- 2
  if(is.null(beta$mu_s2b)) beta$mu_s2b <- 2

  if(is.null(sigma2$a)) sigma2$a <- 0.5
  if(is.null(sigma2$b)) sigma2$b <- 0.5
  if(is.null(sigma2$u)) sigma2$u <- 2

  if(is.null(nu$a)) nu$a <- 1/30
  if(is.null(nu$u)) nu$u <- 2

  if(is.null(dist)) dist <- "student"

  return(list(theta = theta, g = g, tau2 = tau2, beta = beta, sigma2 = sigma2,
              nu = nu, ref = ref, nug_onepar = nug_onepar, alpha = alpha,
              joint_update = joint_update, lambda = lambda, m = m,
              lambda_min = lambda_min, c0 = c0, n0 = n0, dist = dist))
}