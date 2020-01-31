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
#' load("data/ssb.Rdata")
#' source("R/qboot.R")
#' source("R/qbootplot.R")
#' X <- array(ssb, dim = dim(ssb)[c(2,6)])
#' dimnames(X) <- dimnames(ssb)[c(2,6)]
#' X <- t(X)
#' res <- qboot(X, nboot = 100, rerrTarg = 1, aggfun = "max")
#' qbootplot(res, pch = 1)
#' 
qbootplot <- function(obj, pch = 1, col = NULL, ...){
  if(missing(col)) col <- seq(obj$quant)
  
  if(is.null(obj$fit)){
    stop("Nothing to plot. No fitting of model conducted: log(niter)~log(rerr)+quant")
  }
  
  newdat <- expand.grid(
    rerr = seq(min(obj$agg$rerr), max(obj$agg$rerr), length.out = 100),
    quant = unique(obj$agg$quant))
  newdat$niter <- exp(predict(obj$fit, newdata = newdat)) 
  
  plot(niter ~ rerr, obj$agg, log = "xy", pch = pch, col = col, ...)
  legend("topright", legend = paste(obj$quant, "=", round(obj$fitpred$niter)), 
    pch = pch, col = col, ...)

  for(i in seq(levels(newdat$quant))){
    lines(niter ~ rerr, data = newdat, 
      subset = quant == levels(newdat$quant)[i], 
      col = col[i], ...)
  }
  abline(v = obj$rerrTarg, col = 8)
  for(i in seq(levels(newdat$quant))){
    hit <- which(factor(obj$fitpred$quant) == levels(newdat$quant)[i])
    segments(x0 = 1e-6, x1 = obj$rerrTarg, 
      y0 = obj$fitpred$niter[hit], y1 = obj$fitpred$niter[hit], 
      col = col[i], lty = 2)
  }
}

  