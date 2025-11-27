# predict methods for "wrapgp" class

# wrapped gp predict method
#' @title Wrapped Gaussian Process Prediction
#' @description
#' Prediction method for wrapped Gaussian Process model.
#' @param object of class `wrapgp`.
#' @param xx `matrix` of test inputs to predict on. Must have same number of columns as
#' input data used for training.
#' @param lite `boolean` indicating whether to draw from full posterior predictive
#' distribution or sample from individual marginal distributions. The later is
#' typically faster but not as precise.
#' @param verb `boolean` indicating whether to output prediction progress.
#' @export
#' @returns a list containing:
#' \item{mu}{posterior predictive mean estimates.}
#' \item{sd}{posterior predictive standard deviation estimates.}
#' \item{lwr}{posterior predictive lower 2.5% percentile.}
#' \item{upr}{posterior predictive upper 97.5% percentile.}
#' \item{zz}{posterior predictive latent GP samples.}
setMethod("predict", "wrapgp", function(object, xx, lite = T, verb = T){
  # extract mcmc results
  x <- object@x
  dx <- object@dx
  mean_struct <- object@mean_struct
  w <- object@w
  k <- object@k
  kmin <- object@k_min
  z <- object@z
  beta <- object@beta
  mu <- object@mu
  theta <- object@theta
  tau2 <- object@tau2
  g <- object@g
  sigma2 <- object@sigma2
  nug_onepar <- object@nug_onepar
  dist <- object@dist
  nu <- object@nu
  method <- object@method

  n <- nrow(x)
  npred <- nrow(xx)
  nmcmc <- nrow(mu)

  # add intercept term to pred inputs
  xxd <- matrix(rep(1, npred))
  if(mean_struct == "linear") xxd <- as.matrix(cbind(xxd, xx))

  # construct pairwise distance matrices
  dxx <- laGP::distance(xx)
  dxxx <- laGP::distance(xx, x)

  # store pred draws
  mus <- sd2s <- yys <- kks <- zzs <- matrix(nrow = nmcmc, ncol = npred)
  for(t in 1:nmcmc){
    # draw zz
    if(nug_onepar){
      Sx <- tau2[t] * sq_exp(dx, theta[t], 0)
      diag(Sx) <- diag(Sx) + g[t]
      Sxi <- solve(Sx)
      Sxx <- tau2[t] * sq_exp(dxx, theta[t], 0)
      diag(Sxx) <- diag(Sxx) + g[t]
      Sxxx <- tau2[t] * sq_exp(dxxx, theta[t], 0)
    } else{
      Sx <- tau2[t] * sq_exp(dx, theta[t], g[t])
      Sxi <- solve(Sx)
      Sxx <- tau2[t] * sq_exp(dxx, theta[t], g[t])
      Sxxx <- tau2[t] * sq_exp(dxxx, theta[t], 0)
    }
    mu_zz <- drop((xxd %*% beta[t,]) + (Sxxx %*% Sxi %*% (z[t,] - mu[t,])))
    Sigma_zz <- Sxx - (Sxxx %*% Sxi %*% t(Sxxx))
    if(lite){
      zz <- sapply(1:npred, function(i) stats::rnorm(1, mu_zz[i], sqrt(Sigma_zz[i, i])))
    }else{
      zz <- drop(mvtnorm::rmvnorm(1, mu_zz, Sigma_zz))
    }
    zzs[t,] <- zz

    if(method == "jona"){
      mus[t,] <- zz %% (2*pi)
      kks[t,] <- zz %/% (2*pi)
      yys[t,] <- mus[t,]
    }else{
      # draw k
      kk <- eval_tree(xx, w[[t]], kmin)
      kks[t,] <- kk
      
      mus[t,] <- zz - 2*pi*kk
      # draw yy
      if(dist == "normal")  yys[t,] <- (zz - 2*pi*kk + stats::rnorm(npred, sd = sqrt(sigma2[t]))) %% (2*pi)
      if(dist == "student") yys[t,] <- (zz - 2*pi*kk + sqrt(sigma2[t])*stats::rt(npred, nu[t])) %% (2*pi)
    }
    
    if(verb == T & (t %% 100 == 0)) cat(paste0("prediction progress: ", t, "/", nrow(mu), " (", round(100*(t / nrow(mu)), 2), "%) \n"))
  }
  # return mean, quantiles, post draws
  if(dist == "normal"){
    s2 <- mean(sigma2) + sapply(1:npred, function(i) sd(circular(zzs[,i] - 2*pi*kks[,i], units = "radians", modulo = "2pi")))
  }else{
    s2 <- mean(sigma2 * (nu/(nu - 2))) + sapply(1:npred, function(i) sd(circular::circular(zzs[,i] - 2*pi*kks[,i], units = "radians", modulo = "2pi")))
  }
  return(list(yy = yys,
              mu = apply(mus, 2, function(i) mean(circular::circular(i, units = "radians", modulo = "2pi"))),
              s2 = s2,
              lwr = apply(yys, 2, function(i) circular::quantile.circular(circular::circular(i, units = "radians", modulo = "2pi"), 0.025)),
              upr = apply(yys, 2, function(i) circular::quantile.circular(circular::circular(i, units = "radians", modulo = "2pi"), 0.975)),
              zz = zzs,
              kk = kks))
})