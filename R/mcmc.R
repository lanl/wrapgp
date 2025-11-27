# functions called in wrapped gp and deep wrapped gp gibbs samplers

# squared-exponential cov construction
sq_exp <- function(dx, theta, g){
  mat <- exp(-dx / theta)
  diag(mat) <- diag(mat) + g
  return(mat)
}

# mvn log like (z)
llike_z <- function(z, mu, x, dx, theta, tau2, g, nug_onepar = F, profile = F,
                  lambda = 1){
	mu <- c(mu)
	n <- nrow(dx)
	if(nug_onepar){
	  Kn <- sq_exp(dx, theta, g / tau2)
	}else{
	  Kn <- sq_exp(dx, theta, g)
	}
	Kni <- solve(Kn)
	if(profile){
	  xKnix <- solve(t(x) %*% Kni %*% x)
	  bhat <- xKnix %*% (t(x) %*% Kni %*% z)
	  ll <- lambda * (0.5 * msos::logdet(xKnix) - 0.5 * msos::logdet(Kn) -
	    0.5 * (n - 3) * log(0.5 * (t(z - x %*% bhat) %*% Kni %*% (z - x %*% bhat))))
	}else{
	  ll <- lambda * (-0.5 * msos::logdet(tau2 * Kn) - 0.5 * (t(z - mu) %*% (Kni / tau2) %*% (z - mu)))
	}
	return(drop(ll))
}

# normal log like (y)
llike_y <- function(y, z, k, sigma2, dist = "normal", nu = 100, lambda = 1){
  n <- length(y)
  mu <- z - 2*pi*k
  if(dist == "normal"){
    return(lambda * sum(stats::dnorm(y, mu, sqrt(sigma2), log = T)))
  }
  if(dist == "student"){
    return(lambda * sum(stats::dt((y - mu) / sqrt(sigma2), nu, log = T) - log(sqrt(sigma2))))
  }
}

# mvn log like (beta)
llike_beta <- function(beta, dd, theta, g, tau2, mu){
  m <- length(beta)
  Km <- sq_exp(dd, theta, g)
  Kmi <- solve(Km)
  ll <- -0.5 * msos::logdet(tau2 * Km) - 0.5 * (t(beta - mu) %*% (Kmi / tau2) %*% (beta - mu))
  return(drop(ll))
}

# lengthscale metropolis-hastings sampler
# prior = "uniform" -> theta ~ U(a0, b0)
# prior = "gamma" -> theta ~ Gamma(shape = a0, rate = b0)
# other priors unsupported
sample_theta <- function(z, mu, x, dx, theta, tau2, g, u = 1.3, prior = "uniform",
                         a0 = 1e-3, b0 = 1e2, nug_onepar = F, ref = F, lambda = 1){
  # propose new theta
  theta_star <- stats::runif(1, (1/u) * theta, u * theta)

  # construct acceptance threshold
  lr <- llike_z(z, mu, x, dx, theta_star, tau2, g, nug_onepar, ref, lambda) -
    llike_z(z, mu, x, dx, theta, tau2, g, nug_onepar, ref, lambda) +
    log(theta) - log(theta_star)
  if(prior == "uniform"){
    lr <- lr + stats::dunif(theta_star, min = a0, max = b0, log = T) - stats::dunif(theta, min = a0, max = b0, log = T)
  } else if(prior == "gamma"){
    lr <- lr + stats::dgamma(theta_star, shape = a0, rate = b0, log = T) - stats::dgamma(theta, shape = a0, rate = b0, log = T)
  } else{
    stop("specified theta prior currently unsupported")
  }

  # accept/reject proposal
  if(lr >= log(stats::runif(1)) && prior == "uniform" && theta_star >= a0 && theta_star <= b0){
    theta <- theta_star
    accept <- 1
  } else if(lr >= log(stats::runif(1)) && prior == "gamma"){
    theta <- theta_star
    accept <- 1
  } else{
    accept <- 0
  }

  return(list(theta = theta, accept = accept))
}

# scale sampler
# this assumes tau2 ~ IG(a0, b0)
sample_tau2 <- function(z, mu, x, dx, theta, g, tau2, a0, b0, nug_onepar = F, ref = F,
                        lambda = 1){
	if(nug_onepar){
		Kn <- sq_exp(dx, theta, g / tau2)
	}else{
		Kn <- sq_exp(dx, theta, g)
	}
	Kni <- solve(Kn)
	n <- length(z)

	if(ref){
	  bhat <- solve(t(x) %*% Kni %*% x) %*% (t(x) %*% Kni %*% z)
	  samp <- (1/n) * (t(z - x %*% bhat) %*% Kni %*% (z - x %*% bhat))
	}else{
	  an <- (1/2)*(a0 + lambda*n)
	  z_resid <- z - mu
	  bn <- (1/2)*(b0 + lambda*drop(t(z_resid) %*% Kni %*% z_resid))
	  samp <- 1 / stats::rgamma(1, shape = an, rate = bn)
	}

	# Sampling from right-truncated distribution
	#right_constraint <- 10
	#F.right <- pscl::pigamma(right_constraint, an, bn)
	#u <- runif(1, min=0, max=F.right)
	#samp <- qigamma(u, an, bn)
	return(drop(samp))
}

# nugget metropolis-hastings sampler
# currently only supports prior = "gamma" -> g ~ Gamma(shape = a_g, rate = b_g)
# lb specifies a lower bound for the proposed nugget
sample_g <- function(z, mu, x, dx, theta, tau2, g, u = 1.01, prior = "gamma",
                     a_g = 20, b_g = 1/2, nug_onepar = F, lb = 1e-9, ref = F, lambda = 1){
  # propose new g
  g_star <- stats::runif(1, (1/u) * g, u * g)

  # construct acceptance threshold
  lr <- llike_z(z, mu, x, dx, theta, tau2, g_star, nug_onepar, ref, lambda) -
    llike_z(z, mu, x, dx, theta, tau2, g, nug_onepar, ref, lambda) +
    log(g) - log(g_star)
  if(prior == "gamma"){
    lr <- lr + stats::dgamma(g_star, shape = a_g, rate = b_g, log = T) -
      stats::dgamma(g, shape = a_g, rate = b_g, log = T)
  } else{
    stop("specified g prior currently unsupported")
  }

  # accept/reject proposal
  if(lr >= log(stats::runif(1)) && g_star >= lb){
    g <- g_star
    accept <- 1
  } else{
    accept <- 0
  }

  return(list(g = g, accept = accept))
}

# mean coefficient estimates sampler
sample_beta <- function(x, z, Kxi, beta, tau2, mu0, Sigma0, ref = F, alpha = 0, lambda = 1){
  if(ref){
    samp <- solve(t(x) %*% Kxi %*% x) %*% t(x) %*% Kxi %*% z
  }else{
    Ai <- tau2 * solve((lambda*(t(x) %*% Kxi %*% x)) + solve(Sigma0))
    b <- (lambda*(t(x) %*% Kxi %*% z) + solve(Sigma0) %*% mu0) / tau2
    mun <- Ai %*% b
    Sigman <- Ai

    # some slight numerical instabilities can lead to non-symmetric matrices
    if(!isSymmetric(Sigman)){
      Sigman[lower.tri(Sigman)] <- t(Sigman)[lower.tri(Sigman)]
    }
    samp <- mvtnorm::rmvnorm(1, mean = mun + alpha * (beta - mun),
                             sigma = (1 - alpha^2) * Sigman)
  }
  return(drop(samp))
}

# sum-of-squares residuals
ssr <- function(z, mu) return(sum((z - mu)^2))

# proposed ess sampling approach for winding numbers (currently in use)
sample_z <- function(x, y, k, z, dx, mu, theta, tau2, g, sigma2, dist, nu,
                     nug_onepar = F, ref = F, lambda = 1, max_iter = 100){
  # draw u, set ll threshold
  u <- log(stats::runif(1))
  ll_curr <- llike_y(y, z, k, sigma2, dist, nu, lambda)
  ll_thresh <- ll_curr + u

  r <- stats::runif(1, 0, 2*pi)
  r_min <- r - 2*pi
  r_max <- r

  # if reference prior remove scale
  if(ref) tau2 <- 1

  # generate proposal
  if(nug_onepar){
    z0 <- drop(mvtnorm::rmvnorm(1, sigma = tau2 * sq_exp(dx, theta, g / tau2)))
  } else{
    z0 <- drop(mvtnorm::rmvnorm(1, sigma = tau2 * sq_exp(dx, theta, g)))
  }
  z_star <- (z - mu)*cos(r) + z0*sin(r)

  ll_star <- llike_y(y, z_star + mu, k, sigma2, dist, nu, lambda)

  # compare lls, tweak angle r if needed
  z_samps <- z_star
  lls <- ll_star
  rs <- r
  num_iter <- 0
  while(num_iter < max_iter){
    z_samps <- rbind(z_samps, z_star + mu)
    lls <- c(lls, ll_star)
    rs <- c(rs, r)
    if(ll_star > ll_thresh){
      return(list(z = z_star + mu, accept = 1, r = r, z_samps = z_samps, lls = lls,
                  rs = rs, z0 = z0))
    }else{
      if(r < 0){
        r_min <- r
      }else{
        r_max <- r
      }
      # tweak proposal, reevaluate ll
      r <- stats::runif(1, r_min, r_max)
      z_star <- (z - mu)*cos(r) + z0*sin(r)
      ll_star <- llike_y(y, z_star + mu, k, sigma2, dist, nu, lambda)
      num_iter <- num_iter + 1
    }
  }
  cat("failed to accept z_prop \n")
  return(list(z = z, accept = 0, r = 0, z_samps = z_samps, lls = lls,
              rs = rs, z0 = z0))
}

eval_tree <- function(x, w, kmin = 0){
  k <- sapply(x, function(i) sum(i >= w))
  return(k + kmin)
}

shift <- function(z, y, x, w, k, kmin, sigma2, dist, nu, lambda){
  # save copies of w, k
  w_curr <- w
  k_curr <- k
  num_wrap <- length(w_curr)
  if(num_wrap == 0){
    return(list(w_curr, k = k_curr))
  }
  if(num_wrap == 1){
    w <- stats::runif(1, min(x), max(x))
    # propose new k
    k <- eval_tree(x, w, kmin)
    # accept/reject proposed shift
    lr <- llike_y(y, z, k, sigma2, dist, nu, lambda) - 
      llike_y(y, z, k_curr, sigma2, dist, nu, lambda)
    if(lr > log(stats::runif(1))){
      return(list(w = w, k = k))
    }else{
      return(list(w = w_curr, k = k_curr))
    }
  }
  for(i in 1:num_wrap){
    # propose shift for wrapping number
    if(i == 1){
      w[i] <- stats::runif(1, min(x), w[i + 1])
    }else{
      if(i == num_wrap){
        w[i] <- stats::runif(1, w[i - 1], max(x))
      }else{
        w[i] <- stats::runif(1, w[i - 1], w[i + 1])
      }
    }
    # propose new k
    k <- eval_tree(x, w, kmin)
    # accept/reject proposed shift
    lr <- llike_y(y, z, k, sigma2, dist, nu, lambda) - 
      llike_y(y, z, k_curr, sigma2, dist, nu, lambda)
    if(lr > log(stats::runif(1))){
      w_curr <- w
      k_curr <- k
    }
  }
  return(list(w = w_curr, k = k_curr))
}

grow <- function(z, y, x, w, k, kmin, sigma2, dist, nu, lambda){
  # save copies of w, k
  w_curr <- w
  k_curr <- k
  num_wrap <- length(w_curr)
  # propose new wrapping number
  w_new <- stats::runif(1, min(x), max(x))
  w <- sort(c(w, w_new))
  # propose new k
  k <- eval_tree(x, w, kmin)
  # accept/reject proposed shift
  lr <- llike_y(y, z, k, sigma2, dist, nu, lambda) - 
    llike_y(y, z, k_curr, sigma2, dist, nu, lambda) +
    abs(max(x) - min(x)) - log(num_wrap + 1)
  if(lr > log(stats::runif(1))){
    return(list(w = w, k = k))
  }else{
    return(list(w = w_curr, k = k_curr))
  }
}

shrink <- function(z, y, x, w, k, kmin, sigma2, dist, nu, lambda){
  # save copies of w, k
  w_curr <- w
  k_curr <- k
  num_wrap <- length(w_curr)
  # remove wrapping number at random
  w <- w[-sample(1:num_wrap, 1)]
  # propose new k
  k <- eval_tree(x, w, kmin)
  # accept/reject proposed shift
  lr <- llike_y(y, z, k, sigma2, dist, nu, lambda) - 
    llike_y(y, z, k_curr, sigma2, dist, nu, lambda) + 
    log(num_wrap) - abs(max(x) - min(x))
  if(lr > log(stats::runif(1))){
    return(list(w = w, k = k))
  }else{
    return(list(w = w_curr, k = k_curr))
  }
}

sample_k <- function(z, y, k, sigma2, dist, nu, lambda){
  n <- length(y)
  # start at leftmost point
  prop_set <- setdiff(c(k[1] - 1, k[2]), k[1])
  if(length(prop_set) == 1){
    k_star <- prop_set
  }else{
    k_star <- sample(prop_set, 1) # propose new k
  }
  lr <- llike_y(y[1], z[1], k_star, sigma2, dist, nu, lambda) - # eval ratio
    llike_y(y[1], z[1], k[1], sigma2, dist, nu, lambda)
  if(lr > log(stats::runif(1))) k[1] <- k_star # accept/reject proposal

  # iterate through all but last point
  for(i in 2:(n - 1)){
    prop_set <- setdiff(c(k[i - 1], k[i + 1]), k[i])
    if(length(prop_set) == 0) next # if nothing can be proposed, skip
    if(length(prop_set) == 1){
      k_star <- prop_set
    }else{
      k_star <- sample(prop_set, 1) # propose new k
    }
    lr <- llike_y(y[i], z[i], k_star, sigma2, dist, nu, lambda) - # eval ratio
      llike_y(y[i], z[i], k[i], sigma2, dist, nu, lambda)
    if(lr > log(stats::runif(1))) k[i] <- k_star # accept/reject proposal
  }

  # end with rightmost point
  prop_set <- setdiff(c(k[n - 1], k[n] + 1), k[n])
  if(length(prop_set) == 1){
    k_star <- prop_set
  }else{
    k_star <- sample(prop_set, 1) # propose new k
  }
  lr <- llike_y(y[n], z[n], k_star, sigma2, dist, nu, lambda) - # eval ratio
    llike_y(y[n], z[n], k[n], sigma2, dist, nu, lambda)
  if(lr > log(stats::runif(1))) k[n] <- k_star # accept/reject proposal

  return(k)
}

# conditional mean/variance calculations
cond_norm <- function(x, z, ind, Sx, Sxi, mu){
  b <- Sx[ind, -ind]
  g <- Sxi[ind, -ind]
  Sxi <- Sxi[-ind, -ind] + (Sxi[-ind, -ind] %*% b %*% g) / (1 - sum(g * b))
  muc <- mu[ind] + (Sx[ind, -ind, drop = F] %*% Sxi %*% (z[-ind] - mu[-ind]))
  tau2c <- Sx[ind, ind] - (Sx[ind, -ind, drop = F] %*% Sxi %*% Sx[-ind, ind, drop = F])
  return(list(muc = drop(muc), tau2c = drop(tau2c)))
}

# gelfand posterior winding number sampler (currently not in use)
sample_k_jona <- function(x, y, z, dx, mu, theta, g, tau2, lambda, k_min, k_max){
  n <- length(y)
  m <- k_min:k_max
  k <- numeric(n)
  Sx <- tau2 * sq_exp(dx, theta, g)
  Sxi <- solve(Sx)
  for(i in 1:n){
      cond <- cond_norm(x, z, i, Sx, Sxi, mu)
      muc <- cond$muc
      tau2c <- cond$tau2c
      probs <- stats::dnorm(y[i] + 2*pi*m, mean = muc, sd = sqrt(tau2c), log = F)
      probs <- probs / sum(probs)
      k[i] <- sample(m, 1, prob = probs)
  }
  return(k)
}

# sample sigma2
sample_sigma2 <- function(z, y, k, sigma2, a, b, u, dist, nu, lambda){
  if(dist == "normal"){
    n <- length(z)
    an <- (1/2)*(a + n)
    bn <- (1/2)*(b + lambda*drop(t(y + 2*pi*k - z) %*% (y + 2*pi*k - z)))
    return(list(sigma2 = 1 / stats::rgamma(1, an, bn), accept = 1))
  }
  if(dist == "student"){
    # propose new sigma2
    sigma2_star <- stats::runif(1, (1/u) * sigma2, u * sigma2)

    # construct acceptance threshold
    lr <- llike_y(y, z, k, sigma2_star, dist, nu) - llike_y(y, z, k, sigma2, dist, nu) +
      stats::dgamma(sigma2_star, shape = a, rate = b, log = T) -
      stats::dgamma(sigma2, shape = a, rate = b, log = T) +
      log(sigma2) - log(sigma2_star)

    # accept/reject proposal
      if(lr >= log(stats::runif(1))){
        sigma2 <- sigma2_star
        accept <- 1
      }else{
        accept <- 0
      }

      return(list(sigma2 = sigma2, accept = accept))
  }
}

# sample nu
sample_nu <- function(z, y, k, sigma2, nu, a, u, lambda){
  # propose new nu
  nu_star <- stats::runif(1, (1/u) * nu, u * nu)

  # construct acceptance threshold
  lr <- llike_y(y, z, k, sigma2, "student", nu_star, lambda) -
    llike_y(y, z, k, sigma2, "student", nu, lambda) +
    stats::dexp(nu_star, a, log = T) - stats::dexp(nu, a, log = T) +
    log(nu) - log(nu_star)

  # accept/reject proposal
  if(lr >= log(stats::runif(1)) && nu_star >= 3){ # left-truncate at 3
    nu <- nu_star
    accept <- 1
  }else{
    accept <- 0
  }

  return(list(nu = nu, accept = accept))
}

# sample lambda
sample_lambda <- function(z, y, k, sigma2, dist, nu, lambdas, i, p_lambda){
  m <- length(lambdas)
  if(i == 1){
    j <- 2
    qij <- log(1)
    qji <- log(0.5)
  }else{
    if(i == m){
      j <- m - 1
      qij <- log(1)
      qji <- log(0.5)
    }else{
      j <- sample(c(i - 1, i + 1), 1)
      qij <- qji <- log(0.5)
    }
  }

  # accept/reject proposal
  lr <- llike_y(y, z, k, sigma2, dist, nu, lambdas[j]) -
    llike_y(y, z, k, sigma2, dist, nu, lambdas[i]) +
    p_lambda[j] - p_lambda[i] +
    qji - qij

  if(lr > log(stats::runif(1))){
    return(list(i = j, accept = T))
  }else{
    return(list(i = i, accept = F))
  }
}

# sample betas
sample_betas <- function(beta, dd, theta_beta, g_beta, tau2_beta, mu_beta,
                         fits, ind, max_iter = 100){
  m <- length(beta)

  # draw u, set ll threshold
  u <- log(stats::runif(1))
  ll_curr <- 0
  for(i in 1:m){
    x <- fits[[i]]@x
    xd <- cbind(rep(1, nrow(x)), x)
    dx <- fits[[i]]@dx
    z <- fits[[i]]@z[ind,]
    theta <- fits[[i]]@theta[ind]
    tau2 <- fits[[i]]@tau2[ind]
    g <- fits[[i]]@g[ind]
    alpha <- fits[[i]]@beta[ind, 1]
    mu <- drop(xd %*% c(alpha, beta[i]))
    ll_curr <- ll_curr + llike_z(z, mu, x, dx, theta, tau2, g)
  }
  ll_thresh <- ll_curr + u

  r <- stats::runif(1, 0, 2*pi)
  r_min <- r - 2*pi
  r_max <- r

  # generate proposal
  beta0 <- drop(mvtnorm::rmvnorm(1, mean = mu_beta,
                                 sigma = tau2_beta*sq_exp(dd, theta_beta, g_beta)))
  beta_star <- beta*cos(r) + beta0*sin(r)

  ll_star <- 0
  for(i in 1:m){
    x <- fits[[i]]@x
    xd <- cbind(rep(1, nrow(x)), x)
    dx <- fits[[i]]@dx
    z <- fits[[i]]@z[ind,]
    theta <- fits[[i]]@theta[ind]
    tau2 <- fits[[i]]@tau2[ind]
    g <- fits[[i]]@g[ind]
    alpha <- fits[[i]]@beta[ind, 1]
    mu_star <- drop(xd %*% c(alpha, beta_star[i]))
    ll_star <- ll_star + llike_z(z, mu_star, x, dx, theta, tau2, g)
  }

  # compare lls, tweak angle r if needed
  beta_samps <- beta_star
  lls <- ll_star
  rs <- r
  num_iter <- 0
  while(num_iter < max_iter){
    beta_samps <- rbind(beta_samps, beta_star)
    lls <- c(lls, ll_star)
    rs <- c(rs, r)
    if(ll_star > ll_thresh){
      return(list(beta = beta_star, accept = 1, r = r, beta_samps = beta_samps,
                  lls = lls, rs = rs, beta0 = beta0))
    }else{
      if(r < 0){
        r_min <- r
      }else{
        r_max <- r
      }
      # tweak proposal, reevaluate ll
      r <- stats::runif(1, r_min, r_max)
      beta_star <- beta*cos(r) + beta0*sin(r)
      ll_star <- 0
      for(i in 1:m){
        x <- fits[[i]]@x
        xd <- cbind(rep(1, nrow(x)), x)
        dx <- fits[[i]]@dx
        z <- fits[[i]]@z[ind,]
        theta <- fits[[i]]@theta[ind]
        tau2 <- fits[[i]]@tau2[ind]
        g <- fits[[i]]@g[ind]
        alpha <- fits[[i]]@beta[ind, 1]
        mu_star <- drop(xd %*% c(alpha, beta_star[i]))
        ll_star <- ll_star + llike_z(z, mu_star, x, dx, theta, tau2, g)
      }
      num_iter <- num_iter + 1
    }
  }
  cat("failed to accept beta_star \n")
  return(list(beta = beta, accept = 0, r = 0, beta_samps = beta_samps,
              lls = lls, rs = rs, beta0 = beta0))
}

# sample theta_beta
sample_theta_beta <- function(beta, dd, theta, g, tau2, mu, a0 = 1, b0 = 2.6, u = 2){
  # propose new theta
  theta_star <- stats::runif(1, (1/u) * theta, u * theta)

  # construct acceptance threshold
  lr <- llike_beta(beta, dd, theta_star, g, tau2, mu) -
    llike_beta(beta, dd, theta, g, tau2, mu) +
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

# sample g_beta
sample_g_beta <- function(beta, dd, theta, g, tau2, mu, a0 = 1, b0 = 2.6, u = 2){
  # propose new g
  g_star <- stats::runif(1, (1/u) * g, u * g)

  # construct acceptance threshold
  lr <- llike_beta(beta, dd, theta, g_star, tau2, mu) -
    llike_beta(beta, dd, theta, g, tau2, mu) +
    stats::dgamma(g_star, a0, b0, log = T) -
    stats::dgamma(g, a0, b0, log = T) +
    log(g) - log(g_star)

  # accept/reject proposal
  if(lr >= log(stats::runif(1))){
    g <- g_star
    accept <- 1
  }else{
    accept <- 0
  }

  return(list(g = g, accept = accept))
}

# sample tau2_beta
sample_tau2_beta <- function(beta, dd, theta, g, mu, a0 = 2, b0 = 1){
  m <- length(beta)
  Km <- sq_exp(dd, theta, g)
  Kmi <- solve(Km)
  am <- (1/2)*(a0 + m)
  bm <- (1/2)*(b0 + drop(t(beta - mu) %*% Kmi %*% (beta - mu)))
  samp <- 1 / stats::rgamma(1, shape = am, rate = bm)

  return(drop(samp))
}

# sample distance
sample_dist <- function(beta, d, dd, theta, g, tau2, beta_beta, mu, ind, u = 0.1){
  m <- length(beta)

  # save current distance
  d_curr <- d[ind]

  # propose new dist
  d_star <- stats::rnorm(1, d_curr, u)

  # reject if outside feasible bounds
  if(d_star < 0 || d_star > 1) return(list(d = d_curr, accept = 0))

  # construct new covariance, mean
  d[ind] <- d_star
  dd_star <- laGP::distance(d)
  mu_star <- drop(cbind(rep(1, m), d) %*% beta_beta)

  # construct acceptance threshold
  lr <- llike_beta(beta, dd_star, theta, g, tau2, mu_star) -
    llike_beta(beta, dd, theta, g, tau2, mu)

  # accept/reject proposal
  if(lr >= log(stats::runif(1))){
    d_curr <- d_star
    accept <- 1
  }else{
    accept <- 0
  }
  return(list(d = d_curr, accept = accept))
}

# mean coefficient estimate for discrepancy gp
sample_beta_beta <- function(beta, d, dd, theta, g, tau2,
                             mu0 = c(10, 15), Sigma0 = diag(2)){
  Kdi <- solve(sq_exp(dd, theta, g))
  Ai <- tau2 * solve((t(d) %*% Kdi %*% d) + solve(Sigma0))
  b <- (t(d) %*% Kdi %*% beta + solve(Sigma0) %*% mu0) / tau2
  mun <- Ai %*% b
  Sigman <- Ai

  # some slight numerical instabilities can lead to non-symmetric matrices
  if(!isSymmetric(Sigman)){
    Sigman[lower.tri(Sigman)] <- t(Sigman)[lower.tri(Sigman)]
  }
  samp <- mvtnorm::rmvnorm(1, mean = mun, sigma = Sigman)
  return(drop(samp))
}

# sample intercept
sample_alpha <- function(z, x, dx, beta, theta, tau2, g, mu0 = 0, s20 = 1){
  n <- length(z)
  Cxi <- solve(sq_exp(dx, theta, g))
  a <- (1/s20) + (1/tau2)*drop(t(rep(1, n)) %*% Cxi %*% rep(1, n))
  b <- (mu0/s20) + (1/tau2)*drop(t(z - beta*x) %*% Cxi %*% rep(1, n))
  return(stats::rnorm(1, b/a, 1/a))
}
