library(wrapgp)
library(readr)
library(tidyverse)
library(janitor)
library(circular)
set.seed(1)

# load in test data
rfid <- read.csv("data/rfid.csv")

# get unique distances
dists <- unique(rfid$dist)
m <- length(dists)

# standardize distances for training
dist_range <- c(0, 375) # approximate length of lab, in inches
rfid <- rfid %>% mutate(dist = (dist - dist_range[1]) / diff(dist_range))
dist_standard <- (dists - dist_range[1]) / diff(dist_range)

# standardize frequencies for training
freq_range <- c(902.75, 927.25) # frequency range
rfid <- rfid %>% mutate(freq = (freq - freq_range[1]) / diff(freq_range))

# separate data into list of training sets
X <- Y <- list()
for(i in 1:m){
  # grab reads for given distance
  df <- rfid %>% filter(dist == dist_standard[i])
  # randomly sample at most 300 observations for training
  if(nrow(df) > 300){
    df <- df[sort(sample(1:nrow(df), 300, replace = F)),]
  }
  # save training X, Y
  X[[i]] <- matrix(df$freq)
  Y[[i]] <- df$phase
}

# set fixed discrepancy param
mu <- 1

# initialize other discrepancy params
beta_init <- sapply(1:m, function(j) slope_est(X[[j]], Y[[j]], 25))
dat <- data.frame(d = dist_standard, delta = log(beta_init / mu))
lm_fit <- lm(delta ~ d, data = dat)
delta_init <- predict(lm_fit)

# fit hierarchical wrapped gp
nmcmc <- 10000
set.seed(1)
fit <- wrapgp_hier(X, Y, dist_standard, beta = beta_init, delta = delta_init, 
                   mu_delta = mu, nmcmc = nmcmc)

# predict delta on fine grid
inds   <- seq(nmcmc/2, nmcmc, 10)
deltas <- fit$delta[inds,]
thetas <- fit$theta[inds]
tau2s  <- fit$tau2[inds]
betas  <- fit$beta[inds,]

dd <- laGP::distance(dist_standard)

npred <- 500
dg <- seq(0, 1, length.out = npred)
ddg <- laGP::distance(dg)

dgd <- laGP::distance(dist_standard, dg)

delta_draws <- matrix(nrow = nrow(deltas), ncol = npred)
mu_draws <- matrix(nrow = nrow(deltas), ncol = npred)
s2_draws <- matrix(nrow = nrow(deltas), ncol = npred)
for(t in 1:nrow(deltas)){
  Sd <- tau2s[t]*wrapgp:::sq_exp(dd, thetas[t], 1e-4)
  Sdi <- solve(Sd)
  Sdg <- tau2s[t]*wrapgp:::sq_exp(ddg, thetas[t], 1e-4)
  Sdgd <- tau2s[t]*wrapgp:::sq_exp(dgd, thetas[t], 1e-6)
  mu <- drop(t(Sdgd) %*% Sdi %*% deltas[t,])
  Sigma <- Sdg - (t(Sdgd) %*% Sdi %*% Sdgd)
  Sigma <- round(Sigma, 3)
  delta_draws[t,] <- drop(mvtnorm::rmvnorm(1, mu, Sigma))
}

ms <- apply(exp(delta_draws), 2, mean)
lwr <- apply(exp(delta_draws), 2, function(i) quantile(i, 0.025))
upr <- apply(exp(delta_draws), 2, function(i) quantile(i, 0.975))

pdf("rfid_slope_est.pdf", width = 5, height = 5)
plot(dg*diff(dist_range) + dist_range[1], ms, type = "l", col = "black", 
     lwd = 2, xlab = "Distance (in)", ylab = "Slope Estimate", ylim = c(10, 26),
     main = "WGP (Ours)", cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
polygon(c(dg*diff(dist_range) + dist_range[1], rev(dg*diff(dist_range) + dist_range[1])), 
        c(lwr, rev(upr)), col = alpha("grey", 0.5), border = F)
points(dists, apply(betas, 2, mean), col = "red", pch = 19, xlab = "Distance (in)", 
       ylab = "Slope Estimate", ylim = c(0, 25.5))
arrows(dists, apply(betas, 2, function(i) quantile(i, 0.025)),
       dists, apply(betas, 2, function(i) quantile(i, 0.975)),
       col = "red", angle = 90, length = 0.1)
arrows(dists, apply(betas, 2, function(i) quantile(i, 0.975)),
       dists, apply(betas, 2, function(i) quantile(i, 0.025)),
       col = "red", angle = 90, length = 0.1)
legend("bottomright", 
       legend = c(expression(e^delta), expression(beta)),
       col = c("black", "black"), pch = c(NA, 19), lwd = c(2, NA), cex = 1.3)
dev.off()

XX <- as.matrix(seq(0, 1, length.out = 500))

for(i in c(1, 10, 15)){
  pred <- predict(trim(wrap_fit$wgps[[i]], nmcmc / 2, 10), XX)
  mu <- pred$mu
  s2 <- pred$s2
  lwr <- (mu - 1.96*sqrt(s2)) %% (2*pi)
  upr <- (mu + 1.96*sqrt(s2)) %% (2*pi)
  
  lwr2 <- upr2 <- numeric(nrow(XX))
  inds1 <- which(upr < mu)
  upr2[inds1] <- upr[inds1]
  lwr2[inds1] <- 0
  upr[inds1] <- 2*pi
  inds2 <- which(lwr > mu)
  lwr2[inds2] <- lwr[inds2]
  upr2[inds2] <- 2*pi
  lwr[inds2] <- 0
  
  pdf(paste0("phase_fit", i, ".pdf"), width = 5, height = 5)
  plot(X[[i]]*diff(freq_range) + freq_range[1], Y[[i]], pch = 20, 
       xlab = "Frequency (MHz)", main = paste(round(dists[i]), "inches"), ylab = "Phase Angle (radians)", 
       cex.main = 1.7, cex.lab = 1.5, cex.axis = 1.3)
  lines(XX*diff(freq_range) + freq_range[1], mu, col = "red", lwd = 2)
  polygon(c(XX*diff(freq_range) + freq_range[1], rev(XX)*diff(freq_range) + freq_range[1]), 
          c(lwr, rev(upr)), col = alpha("red", 0.3), border = F)
  polygon(c(XX*diff(freq_range) + freq_range[1], rev(XX*diff(freq_range) + freq_range[1])), 
          c(lwr2, rev(upr2)), col = alpha("red", 0.3), border = F)
  dev.off()
}