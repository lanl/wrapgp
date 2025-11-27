library(wrapgp)

methods <- c("wgp_ours", "wgp_jonas", "mallasto")

ns <- seq(100, 250, 50)

nreps <- 10

results <- data.frame(method = "gp", n = -1, seed = -1, rmse = -1, crps = -1)
for(method in methods){
  for(n in ns){
    for(seed in 1:nreps){
      pred_name <- paste0("pred/", method, "_n", n, "_seed", seed, ".csv")
      draws_name <- paste0("draws/", method, "_n", n, "_seed", seed, ".csv")
      if(file.exists(pred_name) & file.exists(draws_name)){
        pred <- read.csv(pred_name)
        mu <- pred$m
        s2 <- pred$s2
        YY <- pred$Y
        rmse <- wrap_rmse(YY, mu)
        
        draws <- read.csv(draws_name)
        crps <- wrap_crps(YY, draws)
        
        results <- rbind(results, c(method, n, seed, rmse, crps))
      }
    }
  }
}

write.csv(results[-1,], "results.csv", row.names = F)
