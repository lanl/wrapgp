library(lhs)

holsetal13log <- function(x){
  x <- 2.5*x
  y <- 15*log(1 + x)
  return(y)
}

nreps <- 10
ns <- c(100, 150, 200, 250)
for(seed in 1:nreps){
  for(n in ns){
    set.seed(seed)
    X <- c(0.15*lhs::randomLHS((0.15/0.8)*n, 1), 0.55*lhs::randomLHS((0.55/0.8)*n, 1) + 0.25)
    X <- sort(c(X, 0.1*lhs::randomLHS(n - length(X), 1) + 0.9))
    Z <- holsetal13log(X) + rnorm(n, sd = 0.5)
    Y <- Z %% (2*pi)
    write.csv(data.frame(X = X, Y = Y), paste0("train_n", n, "_seed", seed, ".csv"),
              row.names = F)
  }
}

for(seed in 1:nreps){
  set.seed(seed + 100)
  XX <- sort(lhs::randomLHS(500, 1))
  ZZ <- holsetal13log(XX)
  YY <- ZZ %% (2*pi)
  write.csv(data.frame(X = XX, Y = YY), paste0("test_seed", seed, ".csv"),
            row.names = F)
}

