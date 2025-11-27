# plot methods for "wrapgp" class

# plot method
#' @title Wrapped Gaussian Process Plot
#' @description
#' Plotting method for wrapped Gaussian Process (WGP) model.
#' @param object of class `wrapgp`.
#' @param x `numeric` matrix of inputs.
#' @param y `numeric` vector of outputs.
#' @param type `string` whether plotted in Cartesian ("flat") or polar ("circle") coordinates. 
#' @param pcol `string` color of training points.
#' @param icol `string` color of prediction line and interval.
#' @param cex_lab `numeric` size of axis labels.
#' @param cex_axis `numeric` size of axis tick values.
#' @param fill `boolean` whether interval should be filled in.
#' @param lite `boolean` whether prediction should be pointwise.
#' @export
#' @returns plot of WGP fit
setMethod("plot", "wrapgp", function(object, x, y, type = "flat", pcol = "black", 
                                  icol = "red", cex_lab = 1, 
                                  cex_axis = 0.8, fill = T, lite = T){
  # extract mcmc results
  x <- object@x
  y <- object@y

  # create fine test grid
  xx <- matrix(seq(min(x), max(x), length.out = 500))

  # predict on test set
  pred <- predict(object, xx, verb = F, lite = lite)
  mu <- pred$mu
  lwr <- (mu - 1.96*sqrt(pred$s2)) %% (2*pi)
  upr <- (mu + 1.96*sqrt(pred$s2)) %% (2*pi)

  if(type == "flat"){
    lwr2 <- upr2 <- numeric(nrow(xx))
    inds1 <- which(upr < mu)
    upr2[inds1] <- upr[inds1]
    lwr2[inds1] <- 0
    upr[inds1] <- 2*pi
    inds2 <- which(lwr > mu)
    lwr2[inds2] <- lwr[inds2]
    upr2[inds2] <- 2*pi
    lwr[inds2] <- 0
    
    plot(xx, mu, type = "l", col = icol, lwd = 2, ylim = c(-0.1, 2*pi + 0.1),
         xlab = "X", ylab = "Y", cex.lab = cex_lab, cex.axis = cex_axis)
    if(fill){
      polygon(c(xx, rev(xx)), c(lwr, rev(upr)), col = alpha(icol, 0.3), border = F)
      polygon(c(xx, rev(xx)), c(lwr2, rev(upr2)), col = alpha(icol, 0.3), border = F)
    }else{
      polygon(c(xx, rev(xx)), c(lwr, rev(upr)), col = "white", border = icol)
      polygon(c(xx, rev(xx)), c(lwr2, rev(upr2)), col = "white", border = icol)
    }
    points(x, y, col = pcol, pch = 20)
  }
  if(type == "circle"){
    # scale input
    x <- x - min(x)
    xx <- xx - min(xx)
    
    xn <- x*cos(y)
    yn <- x*sin(y)
    mux <- xx*cos(mu)
    muy <- xx*sin(mu)
    plot(xn, yn, col = pcol, pch = 20, 
         xlim = range(c(xx*cos(upr), xx*cos(lwr), mux)), 
         ylim = range(c(xx*sin(upr), xx*sin(lwr), muy)), 
         xlab = expression(Xcos(Y)), ylab = expression(Xsin(Y)))
    lines(mux, muy, col = icol, lwd = 2)
    lines(xx*cos(lwr), xx*sin(lwr), col = icol, lwd = 2, lty = 3)
    lines(xx*cos(upr), xx*sin(upr), col = icol, lwd = 2, lty = 3)
    points(xn, yn, col = pcol, pch = 20)
  }
})
