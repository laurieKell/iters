

qbootplot <- function(obj, pch = 1, col = NA, ...){
  if(is.na(col)) col <- seq(obj$quant)
  
  newdat <- expand.grid(
    semu = seq(min(obj$agg$semu), max(obj$agg$semu), length.out = 100),
    quant = unique(obj$agg$quant))
  newdat$niter <- exp(predict(obj$fit, newdata = newdat)) 
  
  plot(niter ~ semu, obj$agg, log = "xy", pch = pch, col = col, ...)
  legend("topright", legend = paste(obj$quant, "=", round(obj$nsim$niter)), 
    pch = pch, col = col)

  for(i in seq(levels(newdat$quant))){
    lines(niter ~ semu, data = newdat, 
      subset = quant == levels(newdat$quant)[i], 
      col = col[i])
  }
  abline(v = obj$targsemu, col = 8)
  for(i in seq(levels(newdat$quant))){
    hit <- which(factor(obj$nsim$quant) == levels(newdat$quant)[i])
    segments(x0 = 1e-6, x1 = obj$targsemu, 
      y0 = obj$nsim$niter[hit], y1 = obj$nsim$niter[hit], 
      col = col[i], lty = 2)
  }
}

  