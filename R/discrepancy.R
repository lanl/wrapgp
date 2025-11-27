# helper functions called by gibbs_hier_wrapped_gp for modeling distance discrepancy

# mvn log like (delta)
llike_delta <- function(delta, dd, theta, tau2){
  m <- length(delta)
  Km <- sq_exp(dd, theta, 1e-4)
  Kmi <- solve(Km)
  ll <- -0.5 * msos::logdet(tau2 * Km) - 0.5 * (t(delta) %*% (Kmi / tau2) %*% delta)
  return(drop(ll))
}

# sample s2
sample_s2 <- function(beta, mu, delta, a0 = 1, b0 = 1){
  m <- length(beta)
  return(1 / stats::rgamma(1, 0.5*(a0 + m), 0.5*(b0 + sum((beta - mu*exp(delta))^2))))
}

# sample delta discrepancy
sample_delta <- function(delta, beta, s2_beta, dd, theta, tau2, mu, max_iter = 100){
  m <- length(delta)

  # draw u, set ll threshold
  u <- log(stats::runif(1))
  ll_curr <- 0
  for(i in 1:m){
    ll_curr <- ll_curr + stats::dnorm(beta[i], mu*exp(delta[i]), sqrt(s2_beta), log = T)
  }
  ll_thresh <- ll_curr + u

  r <- stats::runif(1, 0, 2*pi)
  r_min <- r - 2*pi
  r_max <- r

  # generate proposal
  delta0 <- drop(mvtnorm::rmvnorm(1, sigma = tau2*sq_exp(dd, theta, 1e-4)))
  delta_star <- delta*cos(r) + delta0*sin(r)
  ll_star <- 0
  for(i in 1:m){
    ll_star <- ll_star + stats::dnorm(beta[i], mu*exp(delta_star[i]), sqrt(s2_beta), log = T)
  }

  # compare lls, tweak angle r if needed
  delta_samps <- delta_star
  lls <- ll_star
  rs <- r
  num_iter <- 0
  while(num_iter < max_iter){
    delta_samps <- rbind(delta_samps, delta_star)
    lls <- c(lls, ll_star)
    rs <- c(rs, r)
    if(ll_star > ll_thresh){
      return(list(delta = delta_star, accept = 1, r = r, delta_samps = delta_samps,
                  lls = lls, rs = rs, delta0 = delta0))
    }else{
      if(r < 0){
        r_min <- r
      }else{
        r_max <- r
      }
      # tweak proposal, reevaluate ll
      r <- stats::runif(1, r_min, r_max)
      delta_star <- delta*cos(r) + delta0*sin(r)
      ll_star <- 0
      for(i in 1:m){
        ll_star <- ll_star + stats::dnorm(beta[i], mu*exp(delta_star[i]), sqrt(s2_beta), log = T)
      }
      num_iter <- num_iter + 1
    }
  }
  cat("failed to accept delta_star \n")
  return(list(delta = delta, accept = 0, r = 0, delta_samps = delta_samps,
              lls = lls, rs = rs, delta0 = delta0))
}

# sample delta mean
sample_mu_delta <- function(beta, delta, dd, theta, tau2, s2, mu0 = 1, s20 = 1){
  a <- (1/s20) + (1/s2)*sum((exp(delta))^2)
  b <- (mu0/s20) + (1/s2)*sum(beta*exp(delta))
  return(stats::rnorm(1, mean = b / a, sd = 1 / sqrt(a)))
}

# sample theta_delta
sample_theta_delta <- function(delta, dd, theta, tau2, a0 = 1, b0 = 2.6, u = 2){
  # propose new theta
  theta_star <- stats::runif(1, (1/u) * theta, u * theta)

  # construct acceptance threshold
  lr <- llike_delta(delta, dd, theta_star, tau2) -
    llike_delta(delta, dd, theta, tau2) +
    stats::dgamma(theta_star, a0, b0, log = T) -
    stats::dgamma(theta, a0, b0, log = T) +
    log(theta) - log(theta_star)

  # accept/reject proposal
  if(lr >= log(stats::runif(1))){
    theta <- theta_star
    accept <- 1
  }else{
    accept <- 0
  }

  return(list(theta = theta, accept = accept))
}

# sample tau2_delta
sample_tau2_delta <- function(delta, dd, theta, a0 = 1, b0 = 1){
  m <- length(delta)
  Kmi <- solve(sq_exp(dd, theta, 1e-4))
  am <- (1/2)*(a0 + m)
  bm <- (1/2)*(b0 + drop(t(delta) %*% Kmi %*% delta))
  samp <- 1 / stats::rgamma(1, shape = am, rate = bm)

  return(drop(samp))
}

# sample missing distance and delta jointly
sample_dist_delta <- function(beta, delta, delta_curr, d, d_curr, dd, Sdi, theta, tau2, mu, 
                              s2, u = 0.05){
  # propose new dist
  d_star <- stats::rnorm(1, d_curr, u)

  # reject if outside feasible bounds
  if(d_star < 0 || d_star > 1) return(list(d = d_curr, delta = delta_curr, accept = 0))
  
  # propose new delta
  ddd <- laGP::distance(d_star, d)
  Sddd <- tau2*sq_exp(ddd, theta, 1e-4)
  mu_star <- drop(Sddd %*% Sdi %*% delta)
  s2_star <- drop(tau2 - Sddd %*% Sdi %*% t(Sddd))
  delta_star <- mu_star

  ddd <- laGP::distance(d_curr, d)
  Sddd <- tau2*sq_exp(ddd, theta, 1e-4)
  mu_curr <- drop(Sddd %*% Sdi %*% delta)
  s2_curr <- drop(tau2 - Sddd %*% Sdi %*% t(Sddd))
  
  # construct acceptance threshold
  lr <- stats::dnorm(beta, mu*exp(delta_star), sqrt(s2), log = T) - 
    stats::dnorm(beta, mu*exp(delta_curr), sqrt(s2), log = T) +
    stats::dnorm(delta_star, 0, sqrt(tau2), log = T) - 
    stats::dnorm(delta_curr, 0, sqrt(tau2), log = T) +
    stats::dnorm(delta_curr, mu_curr, sqrt(s2_curr), log = T) -
    stats::dnorm(delta_star, mu_star, sqrt(s2_star), log = T)

  # accept/reject proposal
  if(lr >= log(stats::runif(1))){
    return(list(d = d_star, delta = delta_star, accept = T))
  }else{
    return(list(d = d_curr, delta = delta_curr, accept = F))
  }
}
