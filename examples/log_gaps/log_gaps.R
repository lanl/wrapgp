library(wrapgp)
nmcmc <- 50000
burn <- nmcmc / 2
thin <- 10

seed <- 1
method <- 1
n <- 250
args <- commandArgs(T)
if(length(args) > 0){
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
}
cat("seed is", seed, "\n")
cat("method is", method, "\n")
cat("n is", n, "\n")
set.seed(seed)

train_name <- paste0("data/train_n", n, "_seed", seed, ".csv")
test_name <- paste0("data/test_seed", seed, ".csv")
train <- read.csv(train_name)
test <- read.csv(test_name)
X <- as.matrix(train$X)
Y <- train$Y
XX <- as.matrix(test$X)
YY <- test$Y

if(!file.exists("pred")) dir.create("./pred")
if(!file.exists("draws")) dir.create("./draws")

if(method == 1){  # our wrapped gp
  beta_init <- slope_est(X, Y, 15)
  settings <- mcmc_settings(X, theta = list(a = 2, b = 1, u = 1.5, prior = "gamma"),
                            beta = list(mean_struct = "linear", mu0 = c(0, beta_init),
                                        Sigma0 = diag(c(10, 10))),
                            m = 3, lambda_min = 0.5, c0 = 100,
                            dist = "student",
                            sigma2 = list(a = 2, b = 1))
  # fit wrapped model
  fit <- fit_wrapped_gp(X, Y, g = 1e-4, kmin = -2, theta = NULL, tau2 = NULL, beta = NULL, sigma2 = NULL,
                        nmcmc = nmcmc, temp = F, settings = settings)
  
  name <- "wgp_ours"
}

if(method == 2){ #jona-lassino's wrapped gp
  beta_init <- slope_est(X, Y, 15)
  settings <- mcmc_settings(X, theta = list(a = 2, b = 1, u = 1.5, prior = "gamma"),
                            beta = list(mean_struct = "linear", mu0 = c(5, beta_init),
                                        Sigma0 = diag(c(10, 10))),
                            m = 3, lambda_min = 0.5, c0 = 100,
                            dist = "normal",
                            sigma2 = list(a = 2, b = 1))
  # fit wrapped model
  fit <- fit_wrapped_gp(X, Y, g = NULL, theta = NULL, tau2 = NULL, beta = NULL, sigma2 = NULL,
                        nmcmc = nmcmc, temp = F, method = "jonas", k_min = -5, k_max = 5, settings = settings)
  
  name <- "wgp_jonas"
}

# trim mcmc chains
fit <- trim(fit, burn, thin)

# predict on test inputs
pred <- predict(fit, XX)

# save predictions and draws
fname <- paste0(name, "_n", n, "_seed", seed, ".csv")
write.csv(data.frame(X = XX, Y = YY, m = pred$mu, s2 = pred$s2), paste0("pred/", fname), row.names = F)
write.csv(data.frame(pred$yy), paste0("draws/", fname), row.names = F)

