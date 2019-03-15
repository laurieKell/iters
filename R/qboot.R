#' Calculate MSE quantile estimation error via bootstrapping and predict
#' number of iterations needed for a given precision.
#' 
#' @description Use of relative standard error ('semu')
#'
#' @param X matrix. Statistic of interest across years (rows) and 
#'   iterations (columns)
#' @param nboot numeric. Number of bootstraps. 
#' @param quant vector. Quantile probabilities 
#'   (Default: \code{quant = c(0.95,0.75,0.5,0.25,0.05)}).
#' @param niter vector. Cumulative number of iterations to perform bootstrapped
#'   estimates (Default: \code{niter = exp(seq(log(ncol(X))*0.5, 
#'   log(ncol(X)), length.out = 20))}).
#' @param estfun function. Function name to be used in estimating 
#'   bootstrapped quantile (Default: \code{estfun = "mean"}). 
#' @param aggfun string. Name of function to apply to quantile relative error by year
#'   for use in log-log model (Default: \code{aggfun = "max"}). Use of "max" is 
#'   consistent with 'Prob3', which is the maximum annual probability of SSB 
#'   dropping below Blim. 'Prob2' would be the average annual probability 
#'   of SSB dropping below Blim (i.e. \code{aggfun = "mean"}).    
#' @param verbose logical. Should progress be printed 
#'   (Default: code\{verbose = TRUE})
#' @param rseTarget numeric. Desired relative standard error (as a percent) 
#'   for which to predict the required number of iterations.
#' @param ci numeric. Confidence interval percent to use for  
#'   lower (\code{cilow}) and upper (\code{ciup}) limits of bootrapped 
#'   quantiles. 
#'
#' @return list. Contains ...
#' @export
#'
#' @examples
#' load("data/SSB.Rdata")
#' source("R/qboot.R")
#' source("R/qbootplot.R")
#' X <- array(ssb, dim = dim(ssb)[c(2,6)])
#' dimnames(X) <- dimnames(ssb)[c(2,6)]
#' res <- qboot(X, nboot = 10, rseTarget = 1)
#' qbootplot(res)
#' 
#' 
#' 
qboot <- function(X, nboot = 20, quant = c(0.95,0.75,0.5,0.25,0.05), 
  niter = exp(seq(log(ncol(X))*0.5, log(ncol(X)), length.out = 20)),
  estfun = "mean", aggfun = "max", rseTarget = 1, ci = 0.95, verbose = TRUE
){

  # add some tests of X dims
  
  # empty matrix for bootstrapped quantiles
  qn <- array(NaN, dim=c(nrow(X), length(quant), nboot),
    dimnames = list(year=seq(nrow(X)), quant=quant, boot=seq(nboot)))
  # empty matrix for std.err and mean results
  est <- se <- cilow <- ciup <- array(NaN, dim=c(nrow(X), length(quant), length(niter)),
    dimnames = list(year=seq(nrow(X)), quant=quant, niter=niter))
  for(i in seq(niter)){
    for(n in seq(nboot)){
      if(nrow(X)==1){
        qn[,,n] <- quantile(X[,sample(ncol(X), niter[i], replace = TRUE)], prob=quant)
      } else {
        qn[,,n] <- t(apply(X[,sample(ncol(X), niter[i], replace = TRUE)], 
          MARGIN = 1, FUN = quantile, prob=quant))
      }
    }
    est[,,i] <- apply(qn, MARGIN = 1:2, FUN = estfun)
    se[,,i] <- apply(qn, MARGIN = 1:2, FUN = sd) # se of estimate
    cilow[,,i] <- apply(qn, MARGIN = 1:2, FUN = quantile, prob = (1-ci)/2 )
    ciup[,,i] <- apply(qn, MARGIN = 1:2, FUN = quantile, prob = (ci+(1-ci)/2) )

    if(verbose) print(paste0("niter = ", i, " of ", length(niter), " is finished"))
  }
  
  rse <- se/est*100 # relative standard error 
  
  # fit log-log linear regression; lm(log(niter) ~ log(semu) + quant)
  df <- expand.grid(dimnames(rse))
  df$rse <- c(rse)
  agg <- aggregate(rse ~ quant + niter, data = df, FUN = aggfun)
  agg$niter <- as.numeric(as.character(agg$niter))
  fit <- lm(log(niter) ~ log(rse) + quant, agg)
  newdat <- data.frame(rse = rseTarget, quant = factor(quant))
  newdat$niter <- exp(predict(fit, newdata = newdat))
  
  # return results
  ret <- list(est = est, se = se, cilow = cilow, ciup = ciup, rse = rse,
    fit = fit, nsim = newdat,
    nboot = nboot, quant = quant, niter = niter, agg = agg,
    aggfun = aggfun, estfun = estfun, ci = ci, rseTarget = rseTarget)
  class(ret) <- "qboot"
  return(ret)
}
