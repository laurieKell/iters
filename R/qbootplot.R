#' Plot results of qboot
#'
#' @param obj 
#' @param pch 
#' @param col 
#' @param ... 
#'
#' @return plot of qboot results
#' @export
#'
#' @examples
#' load("data/SSB.Rdata")
#' source("R/qboot.R")
#' source("R/qbootplot.R")
#' X <- array(ssb, dim = dim(ssb)[c(2,6)])
#' dimnames(X) <- dimnames(ssb)[c(2,6)]
#' res <- qboot(X, nboot = 20, rseTarget = 1, aggfun = "max")
#' qbootplot(res, pch = 1)
#' 
qbootplot <- function(obj, pch = 1, col = NULL, ...){
  if(missing(col)) col <- seq(obj$quant)
  
  newdat <- expand.grid(
    rse = seq(min(obj$agg$rse), max(obj$agg$rse), length.out = 100),
    quant = unique(obj$agg$quant))
  newdat$niter <- exp(predict(obj$fit, newdata = newdat)) 
  
  plot(niter ~ rse, obj$agg, log = "xy", pch = pch, col = col, ...)
  legend("topright", legend = paste(obj$quant, "=", round(obj$nsim$niter)), 
    pch = pch, col = col, ...)

  for(i in seq(levels(newdat$quant))){
    lines(niter ~ rse, data = newdat, 
      subset = quant == levels(newdat$quant)[i], 
      col = col[i], ...)
  }
  abline(v = obj$rseTarget, col = 8)
  for(i in seq(levels(newdat$quant))){
    hit <- which(factor(obj$nsim$quant) == levels(newdat$quant)[i])
    segments(x0 = 1e-6, x1 = obj$rseTarget, 
      y0 = obj$nsim$niter[hit], y1 = obj$nsim$niter[hit], 
      col = col[i], lty = 2)
  }
}

  