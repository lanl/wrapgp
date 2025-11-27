library(lhs)

holsetal13log <- function(x){
  x <- 2.5*x
  y <- 15*log(1 + x)
  return(y)
}

nreps <- 10
ns <- c(50, 100, 150, 200)
for(seed in 1:nreps){
  for(n in ns){
    set.seed(seed)
    X <- sort(lhs::randomLHS(n, 1))
    Z <- holsetal13log(X) + rnorm(n, sd = 0.5) #metRology::rt.scaled(n, 5, 0, 0.25, 0)
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

