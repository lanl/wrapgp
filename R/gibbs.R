# gibbs samplers for wrapgp and wrapgp_hier

# define wgp class
#' @title Wrapped Gaussian Process Object
#' @description
#' Object returned after training Wrapped Gaussian Process model. Contains
#' posterior samples from estimated quantities that can be used for prediction and
#' visualization.
#' @slot x WGP inputs.
#' @slot y WGP outputs.
#' @slot mean_struct Mean structure.
#' @slot dx Pairwise-euclidean distance matrix for x.
#' @slot w Posterior wrapping trees.
#' @slot k Posterior winding number draws.
#' @slot k_min Posterior minimum k draws.
#' @slot z Posterior latent GP draws.
#' @slot beta Posterior coefficient draws.
#' @slot mu Posterior latent mean draws.
#' @slot theta Posterior lengthscale draws.
#' @slot tau2 Posterior scale draws.
#' @slot g Posterior nugget draws.
#' @slot nug_onepar Scaled nugget value.
#' @slot r Accepted angle draws from elliptical slice sampling.
#' @slot dist Likelihood type.
#' @slot ll_y Likelihood values.
#' @slot ll_z Pseudo-likelihood values from Gaussian process prior.
#' @slot lambda (Inverse) tempering parameter.
#' @slot sigma2 Posterior variance inflation draws.
#' @slot nu Posterior degrees of freedom draws.
#' @slot method Wrapping number sampling method.
#' @export
setClass("wrapgp",
         slots = c(
           x = "matrix",
           y = "numeric",
           mean_struct = "character",
           dx = "matrix",
           w = "list",
           k = "matrix",
           k_min = "numeric",
           z = "matrix",
           beta = "matrix",
           mu = "matrix",
           theta = "numeric",
           tau2 = "numeric",
           g = "numeric",
           nug_onepar = "logical",
           r = "numeric",
           dist = "character",
           ll_y = "numeric",
           ll_z = "numeric",
           lambda = "numeric",
           sigma2 = "numeric",
           nu = "numeric",
           method = "character"))

# ordinary wrapped gp gibbs
gibbs_wrapgp <- function(x, y, beta, theta, tau2, g, z, w, sigma2, nu,
                         w_start, method, k_min, k_range, nmcmc, temp, 
                         verb, settings){
  # extract settings
  prior_theta <- settings$theta$prior
  a_theta <- settings$theta$a
  b_theta <- settings$theta$b
  u_theta <- settings$theta$u
  prior_g <- settings$g$prior
  a_g <- settings$g$a
  b_g <- settings$g$b
  u_g <- settings$g$u
  prior_tau2 <- settings$tau2$prior
  a_tau <- settings$tau2$a
  b_tau <- settings$tau2$b
  a_sigma <- settings$sigma2$a
  b_sigma <- settings$sigma2$b
  u_sigma <- settings$sigma2$u
  a_nu <- settings$nu$a
  u_nu <- settings$nu$u
  mean_struct <- settings$beta$mean_struct
  mu0 <- settings$beta$mu0
  Sigma0 <- settings$beta$Sigma0
  ref <- settings$ref
  nug_onepar <- settings$nug_onepar
  alpha <- settings$alpha
  joint_update <- settings$joint_update
  dist <- settings$dist

  if(temp){
    m <- settings$m
    lambda_min <- settings$lambda_min
    c0 <- settings$c0
    n0 <- settings$n0
  }else{
    lambda <- settings$lambda
  }

  # get number of rows
  n <- nrow(x)

  # construct design matrix with intercept
  xd <- matrix(rep(1, n))
  if(mean_struct == "linear") xd <- cbind(xd, x)

  # get number of columns (including added intercept)
  d <- ncol(xd)

  # compute pairwise distances for covariance construction
  dx <- laGP::distance(x)

  # initialize beta, mu, theta, tau2, g
  theta_draws <- numeric(nmcmc)
  if(is.null(theta)){
    if(prior_theta == "uniform"){
      theta <- a_theta + 0.5 * (b_theta - a_theta)
    }else{
      theta <- 0.5
    }
    est_theta <- T
  }else{
    est_theta <- F
  }
  theta_draws[1] <- theta
  g_draws <- numeric(nmcmc)
  if(is.null(g)){
    g <- 0.1
    est_g <- T
  }else{
    est_g <- F
  }
  g_draws[1] <- g
  tau2_draws <- numeric(nmcmc)
  if(is.null(tau2)){
    if(prior_tau2 == "uniform"){
      tau2 <- a_tau + 0.5 * (b_tau - a_tau)
    }else{
      tau2 <- 0.1
    }
    est_tau2 <- T
  }else{
    est_tau2 <- F
  }
  tau2_draws[1] <- tau2

  # initialize sigma2
  sigma2_draws <- numeric(nmcmc)
  if(is.null(sigma2)){
    est_sigma2 <- T
    sigma2 <- 1e-1
  }else{
    est_sigma2 <- F
  }
  sigma2_draws[1] <- sigma2

  # initial df (student-t only)
  nu_draws <- numeric(nmcmc)
  if(is.null(nu) & dist == "student"){
    est_nu <- T
    nu <- 30
  }else{
    est_nu <- F
  }
  if(dist == "normal") nu <- 30
  nu_draws[1] <- nu
  
  w_draws <- list()
  k_draws <- matrix(nrow = nmcmc, ncol = n)
  if(is.null(w)){
    est_k <- T
    if(is.null(w_start)){
      w <- c(0.5)
      k <- eval_tree(x, w, k_min)
    }else{
      w <- w_start
      k <- eval_tree(x, w, k_min)
    }
  }else{
    est_k <- F
    k <- eval_tree(x, w, k_min)
  }
  w_draws[[1]] <- w
  k_draws[1,] <- k

  # initialize latent unwrapped GP and corresponding winding numbers
  beta_draws <- matrix(nrow = nmcmc, ncol = d)
  if(is.null(beta)){
    beta <- mu0
    est_beta <- T
  }else{
    est_beta <- F
  }
  beta_draws[1,] <- beta
  mu_draws <- matrix(nrow = nmcmc, ncol = n)
  mu <- drop(xd %*% beta)
  mu_draws[1,] <- mu

  z_draws <- matrix(nrow = nmcmc, ncol = n)
  if(is.null(z)){
    est_z <- T
    if(method == "ess"){
      z <- drop(mvtnorm::rmvnorm(1, mu, tau2 * sq_exp(dx, theta, g)))
    }
    if(method == "jona"){
      z <- y + 2*pi*k
    }
  }else{
    est_z <- F
  }
  z_draws[1,] <- z

  # keep track of log-like, log-post
  ll_y <- numeric(nmcmc)
  ll_z <- numeric(nmcmc)

  ll_y[1] <- llike_y(y, z, k, sigma2, dist, nu)
  ll_z[1] <- llike_z(z, mu, x, dx, theta, tau2, g, nug_onepar = nug_onepar,
                     profile = ref)

  # keep track of lengthscale and/or nugget acceptance rates
  theta_accept <- 0
  tau2_accept <- 0
  g_accept <- 0
  sigma2_accept <- 0
  nu_accept <- 0
  r_final <- rep(NA, nmcmc)

  # save temperatures if temp = T
  if(temp){
    lambdas   <- calc_temp(m, lambda_min)      # construct temperature ladder (lambdas)
    p_lambda  <- log(pseudo_prior(lambdas))    # calculate (log) pseudo-prior for temps

    # initialize lambda (index)
    i <- 1
    # store temperature draws
    lambda_draws  <- numeric(nmcmc)
    lambda_draws[1] <- lambdas[i]
    # keep track of temperature acceptance rate
    lambda_accept <- 0
  }

  # start gibbs
  for(t in 2:nmcmc){
    # set k_min based on smallest z drawn
    if(t == 1000 & method == "ess"){
      k_min <- min(z) %/% (2*pi)
      k <- eval_tree(x, w, k_min)
      cat(paste("minimum wrapping number:", k_min, "\n"))
    }
    
    # set new temp
    if(temp) lambda <- lambdas[i]

    # sample beta
    if(nug_onepar){
      Kxi <- solve(sq_exp(dx, theta, g/tau2))
    }else{
      Kxi <- solve(sq_exp(dx, theta, g))
    }
    if(est_beta){
      beta <- sample_beta(xd, z, Kxi, beta, tau2, mu0, Sigma0, ref, alpha, lambda)
    }
    mu <- xd %*% beta

    # sample theta
    if(est_theta & !joint_update){
      samp <- sample_theta(z, mu, x, dx, theta, tau2, g, u_theta, prior_theta, a_theta, b_theta,
                           nug_onepar, ref, lambda)
      theta <- samp$theta
      theta_accept <- theta_accept + samp$accept
    }

    # sample g
    if(est_g){
      samp <- sample_g(z, mu, x, dx, theta, tau2, g, u_g, prior_g, a_g, b_g,
                       nug_onepar, ref = ref, lambda = lambda)
      g <- samp$g
      g_accept <- g_accept + samp$accept
    }

    # sample tau2
    if(est_tau2 & !joint_update){
      tau2 <- sample_tau2(z, mu, x, dx, theta, g, tau2, a_tau, b_tau, nug_onepar, ref, lambda)
      tau2_accept <- tau2_accept + 1
    }

    # sample z
    if(est_z){
      if(method == "ess"){
        samp <- sample_z(x, y, k, z, dx, mu, theta, tau2, g, sigma2, dist, nu,
                         nug_onepar, ref, lambda)
        z <- samp$z
        r <- samp$r
      }
      if(method == "jona"){
        z <- y + 2*pi*k
        r <- 0
      }
    }else{
      r <- 0
    }

    # sample k
    if(est_k){
      if(method == "ess"){
        # shift
        samp <- shift(z, y, x, w, k, k_min, sigma2, dist, nu, lambda)
        w <- samp$w
        k <- samp$k
        
        # grow
        samp <- grow(z, y, x, w, k, k_min, sigma2, dist, nu, lambda)
        w <- samp$w
        k <- samp$k
        
        # shrink
        samp <- shrink(z, y, x, w, k, k_min, sigma2, dist, nu, lambda)
        w <- samp$w
        k <- samp$k
      }
      if(method == "jona"){
        k <- sample_k_jona(x, y, z, dx, mu, theta, g, tau2, lambda, 
                           k_range[1], k_range[2])
      }
    }

    # sample sigma2
    if(est_sigma2){
      samp <- sample_sigma2(z, y, k, sigma2, a_sigma, b_sigma, u_sigma, dist, nu, lambda)
      sigma2 <- samp$sigma2
      sigma2_accept <- sigma2_accept + samp$accept
    }

    # sample nu
    if(est_nu){
      samp <- sample_nu(z, y, k, sigma2, nu, a_nu, u_nu, lambda)
      nu <- samp$nu
      nu_accept <- nu_accept + samp$accept
    }

    if(temp){
      # propose new lambda
      samp <- sample_lambda(z, y, k, sigma2, dist, nu, lambdas, i, p_lambda)
      i <- samp$i
      lambda_accept <- lambda_accept + samp$accept

      # update pseudo-prior
      p_lambda[i] <- p_lambda[i] - (c0 / (t + n0))
      p_lambda[-i] <- p_lambda[-i] + (c0 / (m * (t + n0)))
    }


    # store log-likes
    ll_y[t] <- llike_y(y, z, k, sigma2, dist, nu)
    ll_z[t] <- llike_z(z, mu, x, dx, theta, tau2, g, nug_onepar = nug_onepar,
                       profile = ref)

    # save draws
    if(temp) lambda_draws[t] <- lambdas[i]
    r_final[t] <- r
    w_draws[[t]] <- w
    k_draws[t,] <- k
    z_draws[t,] <- z
    beta_draws[t,] <- beta
    mu_draws[t,] <- mu
    theta_draws[t] <- theta
    g_draws[t] <- g
    tau2_draws[t] <- tau2
    sigma2_draws[t] <- sigma2
    nu_draws[t] <- nu

    # print progress
    if(verb & (t %% 1000 == 0)) cat(paste0("mcmc progress: ", t, "/", nmcmc,
                                           " (", round(100 * (t / nmcmc), 2), "%) \n"))
  }

  if(verb & est_theta)                       cat(paste0("Lengthscale acceptance rate: ", theta_accept / nmcmc, " \n"))
  if(verb & est_g)                           cat(paste0("Nugget accceptance rate: ", g_accept / nmcmc, " \n"))
  if(verb & est_sigma2 & dist == "student")  cat(paste0("Error variance acceptance rate: ", sigma2_accept / nmcmc), " \n")
  if(verb & est_nu)                          cat(paste0("Degrees of freedom acceptance rate: ", nu_accept / nmcmc), " \n")
  if(verb & temp)                            cat(paste0("Lambda acceptance rate: ", lambda_accept / nmcmc, " \n"))

  if(temp) inds <- which(lambda_draws == 1)
  else inds <- 1:nmcmc

  # return object of class "wrapgp"
  wrapped_gp <- methods::new("wrapgp", x = x, y = y, mean_struct = mean_struct,
                             dx = dx, w = w_draws[inds], k = k_draws[inds,, drop = F],
                             k_min = k_min, z = z_draws[inds,, drop = F],
                             beta = beta_draws[inds,, drop = F],
                             mu = mu_draws[inds,, drop = F],
                             theta = theta_draws[inds],
                             tau2 = tau2_draws[inds],
                             g = g_draws[inds], nug_onepar = nug_onepar,
                             dist = dist, ll_y = ll_y[inds], ll_z = ll_z[inds],
                             sigma2 = sigma2_draws[inds],
                             nu = nu_draws[inds], method = method)

  return(wrapped_gp)
}

# hierarchical gibbs
gibbs_wrapgp_hier <- function(x, y, d, betas, g, delta, mu_delta, 
                              theta_delta, tau2_delta, s2_beta, d_miss, 
                              delta_miss, method, k_min, k_range, nmcmc, 
                              verb, settings){
  # extract settings
  prior_theta <- settings$theta$prior
  a_theta <- settings$theta$a
  b_theta <- settings$theta$b
  u_theta <- settings$theta$u
  prior_g <- settings$g$prior
  a_g <- settings$g$a
  b_g <- settings$g$b
  u_g <- settings$g$u
  prior_tau2 <- settings$tau2$prior
  a_tau <- settings$tau2$a
  b_tau <- settings$tau2$b
  a_sigma <- settings$sigma2$a
  b_sigma <- settings$sigma2$b
  u_sigma <- settings$sigma2$u
  a_nu <- settings$nu$a
  u_nu <- settings$nu$u
  mu_s2a <- settings$beta$mu_s2a
  mu_s2b <- settings$beta$mu_s2b
  ref <- settings$ref
  nug_onepar <- settings$nug_onepar
  alpha <- settings$alpha
  joint_update <- settings$joint_update
  dist <- settings$dist

  m <- length(x)               # number of tests
  n <- unlist(lapply(x, nrow)) # number of obs for each test

  # construct design matrix with intercept
  dx <- list()
  xd <- list()
  for(j in 1:m){
    dx[[j]] <- laGP::distance(x[[j]])
    xd[[j]] <- cbind(rep(1, n[j]), x[[j]])
  }

  if(is.null(g)){
    g <- 1e-1
    est_g <- T
  }else{
    est_g <- F
  }
  
  # initialize test-level parameters
  if(is.null(betas)){
    betas <- numeric(m)
    init_beta <- T
  }else{
    init_beta <- F
  }
  wgps <- list()
  for(j in 1:m){
    # initialize parameters
    theta_draws <- numeric(nmcmc)
    if(prior_theta == "uniform"){
      theta <- a_theta + 0.5 * (b_theta - a_theta) # so that we don't initialize outside support
    }else{
      theta <- 0.5
    }
    theta_draws[1] <- theta
    
    g_draws <- numeric(nmcmc)
    g_draws[1] <- g
    
    tau2_draws <- numeric(nmcmc)
    if(prior_tau2 == "uniform"){
      tau2 <- a_tau + 0.5 * (b_tau - a_tau) # so that we don't initialize outside support
    }else{
      tau2 <- 0.1
    }
    tau2_draws[1] <- tau2

    sigma2_draws <- numeric(nmcmc)
    sigma2 <- 1e-1
    sigma2_draws[1] <- sigma2
    
    nu_draws <- numeric(nmcmc)
    nu <- 30
    nu_draws[1] <- nu
    
    w_draws <- list()
    k_draws <- matrix(nrow = nmcmc, ncol = n[j])
    w <- c(0.5)
    k <- eval_tree(x[[j]], w, k_min)
    k_draws[1,] <- k
    w_draws[[1]] <- w

    beta_draws <- matrix(nrow = nmcmc, ncol = 2)
    if(init_beta){
      beta <- c(0, slope_est(x[[j]], y[[j]], 25))
      betas[j] <- beta[2]
    }else{
      beta <- c(0, betas[j])
    }
    beta_draws[1,] <- beta
    mu_draws <- matrix(nrow = nmcmc, ncol = n[j])
    mu <- drop(xd[[j]] %*% beta)
    mu_draws[1,] <- mu

    z_draws <- matrix(nrow = nmcmc, ncol = n[j])
    z <- drop(mvtnorm::rmvnorm(1, mu, tau2 * sq_exp(dx[[j]], theta, g)))
    #z <- z - min(z)
    z_draws[1,] <- z

    # keep track of log-like, log-post
    ll_y <- numeric(nmcmc)
    ll_z <- numeric(nmcmc)

    ll_y[1] <- llike_y(y[[j]], z, k, sigma2, dist, nu)
    ll_z[1] <- llike_z(z, mu, x[[j]], dx[[j]], theta, tau2, g,
                       nug_onepar = nug_onepar, profile = ref)

    wgps[[j]] <- methods::new("wrapgp", x = x[[j]], y = y[[j]],
                              mean_struct = "linear", dx = dx[[j]],
                              w = w_draws, k = k_draws, k_min = k_min, z = z_draws, 
                              beta = beta_draws, mu = mu_draws, theta = theta_draws,
                              tau2 = tau2_draws, g = g_draws,
                              nug_onepar = nug_onepar, dist = dist, ll_y = ll_y,
                              ll_z = ll_z, sigma2 = sigma2_draws, nu = nu_draws,
                              method = method)
  }

  # store slope draws at hierarchical level
  betas_draws <- matrix(nrow = nmcmc, ncol = m)
  betas_draws[1,] <- betas

  # initialize hierarchical parameters
  mu_deltas <- numeric(nmcmc)
  if(is.null(mu_delta)){
    mu_delta <- 1
    est_mu_delta <- T
  }else{
    est_mu_delta <- F
  }
  mu_deltas[1] <- mu_delta
  
  # discrepancy terms
  deltas <- matrix(nrow = nmcmc, ncol = length(d))
  if(is.null(delta)) delta <- rep(mu_delta, m)
  deltas[1,] <- delta
  
  # identify missing distances to impute (if any)
  num_miss <- length(x) - length(d)
  d_miss_draws <- delta_miss_draws <- matrix(nrow = nmcmc, ncol = num_miss)
  if(num_miss > 0){
    if(is.null(d_miss)){
      d_miss <- numeric(num_miss)
      for(l in 1:num_miss){
        d_miss[l] <- 0.5
      }
    }
    if(is.null(delta_miss)){
      delta_miss <- numeric(num_miss)
      for(l in 1:num_miss){
        delta_miss[l] <- mean(delta)
      }
    }
    d_miss_draws[1,] <- d_miss
    delta_miss_draws[1,] <- delta_miss
    est_d <- T
    delta_full <- c(delta, delta_miss)
  }else{
    est_d <- F
    delta_full <- delta
  }
  
  # calculate distance mat for observed distances
  dd <- laGP::distance(d)
  
  # common variance term for betas
  s2_betas <- numeric(nmcmc)
  if(is.null(s2_beta)){
    s2_beta <- 0.1
    est_s2_beta <- T
  }else{
    est_s2_beta <- F
  }
  s2_betas[1] <- s2_beta
  
  # discrepancy lengthscale
  theta_deltas <- numeric(nmcmc)
  if(is.null(theta_delta)){
    theta_delta <- 0.1
    est_theta_delta <- T
  }else{
    est_theta_delta <- F
  }
  theta_deltas[1] <- theta_delta
  
  # discrepancy scale
  tau2_deltas <- numeric(nmcmc)
  if(is.null(tau2_delta)){
    tau2_delta <- 1
    est_tau2_delta <- T
  }else{
    est_tau2_delta <- F
  }
  tau2_deltas[1] <- tau2_delta
  
  # keep track of lengthscale and/or nugget acceptance rates
  theta_accept <- 0
  tau2_accept <- 0
  g_accept <- 0
  sigma2_accept <- 0
  nu_accept <- 0
  r_final <- rep(NA, nmcmc)
  theta_delta_accept <- 0
  dist_accept <- 0

  Sdi <- (1/tau2_delta) * solve(sq_exp(dd, theta_delta, 1e-4))
  
  # start gibbs
  for(t in 2:nmcmc){
    for(j in 1:m){
      wgp <- wgps[[j]]
      x <- wgp@x
      xd <- cbind(rep(1, nrow(x)), x)
      y <- wgp@y
      dx <- wgp@dx
      theta <- wgp@theta[t - 1]
      g <- wgp@g[t - 1]
      beta <- wgp@beta[t - 1,]
      tau2 <- wgp@tau2[t - 1]
      sigma2 <- wgp@sigma2[t - 1]
      nu <- wgp@nu[t - 1]
      z <- wgp@z[t - 1,]
      w <- wgp@w[[t - 1]]
      k <- wgp@k[t - 1,]
      k_min <- wgp@k_min
      mu <- wgp@mu[t - 1,]
      
      # set k_min based on smallest z drawn
      if(method == "ess" & t == 1000){
        k_min <- min(z) %/% (2*pi)
        wgps[[j]]@k_min <- k_min
        k <- eval_tree(x, w, k_min)
      }

      # sample alpha/beta
      if(nug_onepar){
        Kxi <- solve(sq_exp(dx, theta, g/tau2))
      }else{
        Kxi <- solve(sq_exp(dx, theta, g))
      }
      beta <- sample_beta(xd, z, Kxi, beta, tau2,
                          mu0 = c(0, mu_delta*exp(delta_full[j])), 
                          Sigma0 = diag(c(10, s2_beta)))
      betas[j] <- beta[2]
      mu <- drop(xd %*% beta)
      
      # sample theta
      samp <- sample_theta(z, mu, x, dx, theta, tau2, g, u_theta, prior_theta, 
                           a_theta, b_theta, nug_onepar, ref, 1)
      theta <- samp$theta
      theta_accept <- theta_accept + samp$accept
      
      # sample g
      if(est_g){
        samp <- sample_g(z, mu, x, dx, theta, tau2, g, u_g, prior_g, a_g, b_g,
                         F, ref, 1)
        g <- samp$g
        g_accept <- g_accept + samp$accept
      }
      
      # sample tau2
      tau2 <- sample_tau2(z, mu, x, dx, theta, g, tau2, a_tau, b_tau, nug_onepar, ref, 1)
      tau2_accept <- tau2_accept + 1

      # sample z
      if(method == "ess"){
        samp <- sample_z(x, y, k, z, dx, mu, theta, tau2, g, sigma2, dist, nu,
                         nug_onepar, ref, 1)
        z <- samp$z
        r <- samp$r
      }
      if(method == "jona"){
        z <- y + 2*pi*k
        r <- 0
      }
      
      # sample k
      if(method == "ess"){
        # shift
        samp <- shift(z, y, x, w, k, k_min, sigma2, dist, nu, 1)
        w <- samp$w
        k <- samp$k
        
        # grow
        samp <- grow(z, y, x, w, k, k_min, sigma2, dist, nu, 1)
        w <- samp$w
        k <- samp$k
        
        # shrink
        samp <- shrink(z, y, x, w, k, k_min, sigma2, dist, nu, 1)
        w <- samp$w
        k <- samp$k
      }
      if(method == "jona"){
        k <- sample_k_jona(x, y, z, dx, mu, theta, g, tau2, 1, 
                           k_range[1], k_range[2])
      }
      
      # sample sigma2
      samp <- sample_sigma2(z, y, k, sigma2, a_sigma, b_sigma, u_sigma, dist, nu, 1)
      sigma2 <- samp$sigma2
      sigma2_accept <- sigma2_accept + samp$accept
      
      # sample nu
      samp <- sample_nu(z, y, k, sigma2, nu, a_nu, u_nu, 1)
      nu <- samp$nu
      nu_accept <- nu_accept + samp$accept

      # store log-likes
      wgps[[j]]@ll_y[t] <- llike_y(y, z, k, sigma2, dist, nu)
      wgps[[j]]@ll_z[t] <- llike_z(z, mu, x, dx, theta, tau2, g,
                                   nug_onepar = nug_onepar, profile = ref)

      # save draws
      wgps[[j]]@r[t] <- r
      wgps[[j]]@w[[t]] <- w
      wgps[[j]]@k[t,] <- k
      wgps[[j]]@z[t,] <- z
      wgps[[j]]@beta[t,] <- beta
      wgps[[j]]@mu[t,] <- mu
      wgps[[j]]@theta[t] <- theta
      wgps[[j]]@g[t] <- g
      wgps[[j]]@tau2[t] <- tau2
      wgps[[j]]@sigma2[t] <- sigma2
      wgps[[j]]@nu[t] <- nu
    }

    # sample hierarchical params
    # sample discrepancy for only observed distances
    samp <- sample_delta(delta, betas[1:length(d)], s2_beta, dd, 
                         theta_delta, tau2_delta, mu_delta)
    delta <- samp$delta
    
    # sample mu
    if(est_mu_delta){
      mu_delta <- sample_mu_delta(betas, delta, dd, theta_delta,
                                  tau2_delta, s2_beta)
    }
    
    # sample theta
    if(est_theta_delta){
      samp <- sample_theta_delta(delta, dd, theta_delta, tau2_delta)
      theta_delta <- samp$theta
      theta_delta_accept <- theta_delta_accept + samp$accept
    }
    
    # sample tau2
    if(est_tau2_delta){
      tau2_delta <- sample_tau2_delta(delta, dd, theta_delta)
    }
    
    # sample s2_beta
    if(est_s2_beta){
      s2_beta <- sample_s2(betas[1:length(d)], mu_delta, delta)
    }
    
    if(est_theta_delta | est_tau2_delta){
      Sdi <- (1/tau2_delta) * solve(sq_exp(dd, theta_delta, 1e-2))
    }
    
    # sample missing dists
    if(est_d){
      for(l in 1:num_miss){
        samp <- sample_dist_delta(betas[length(d) + l], delta, delta_miss[l], 
                                  d, d_miss[l], dd, Sdi, theta_delta, tau2_delta, 
                                  mu_delta, s2_beta)
        dist_accept <- dist_accept + samp$accept
        d_miss[l] <- samp$d
        delta_miss[l] <- samp$delta
        delta_full <- c(delta, delta_miss)
      }
    }else{
      delta_full <- delta
    }
    
    # store hierarchical draws
    betas_draws[t,] <- betas
    deltas[t,] <- delta
    mu_deltas[t] <- mu_delta
    theta_deltas[t] <- theta_delta
    tau2_deltas[t] <- tau2_delta
    s2_betas[t] <- s2_beta
    if(est_d){
      d_miss_draws[t,] <- d_miss
      delta_miss_draws[t,] <- delta_miss
    }

    # print progress
    if(verb & (t %% 100 == 0)) cat(paste0("mcmc progress: ", t, "/", nmcmc,
                                           " (", round(100 * (t / nmcmc), 2), "%) \n"))
  }

  if(verb)                      cat(paste0("Lengthscale acceptance rate: ", theta_accept / (nmcmc*m), " \n"))
  if(verb & dist == "student")  cat(paste0("Error variance acceptance rate: ", sigma2_accept / (nmcmc*m), " \n"))
  if(verb)                      cat(paste0("Degrees of freedom acceptance rate: ", nu_accept / (nmcmc*m), " \n"))
  if(verb & est_g)              cat(paste0("Nugget acceptance rate: ", g_accept / (nmcmc*m), " \n"))
  if(verb & est_theta_delta)    cat(paste0("Delta lengthscale acceptance rate: ", theta_delta_accept / nmcmc), " \n")
  if(verb & est_d)              cat(paste0("Distance acceptance rate: ", dist_accept / (nmcmc*(m - length(d)))), " \n")
  
  return(list(wgps = wgps, beta = betas_draws, delta = deltas,
              mu = mu_deltas, theta = theta_deltas, tau2 = tau2_deltas,
              s2_beta = s2_betas, d_miss = d_miss_draws, 
              delta_miss = delta_miss_draws))
}