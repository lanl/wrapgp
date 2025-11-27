# calculate temperatures
calc_temp <- function(m, t_min, spacing = "geometric"){
  dt <- t_min^(1 / (1 - m)) - 1
  t <- numeric(m)
  t[1] <- 1
  if(spacing == "geometric"){
    for(i in 2:m){
      t[i] <- (1 + dt)^(1 - i)
    }
  }
  return(t)
}

# pseudo-prior specification
pseudo_prior <- function(t){
  return(rep(1, length(t)) / length(t))
}
