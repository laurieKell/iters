#' Calculate MSE quantile estimation error via bootstrapping and predict
#' number of iterations needed for a given precision.
#' 
#' @description Use of relative standard error ('semu')
#'
#' @param X matrix. Statistic of interest across years (columns) and 
#'   iterations (rows)
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
#' @param rerrTarg numeric. Desired relative standard error (as a percent) 
#'   for which to predict the required number of iterations.
#' @param ci numeric. Confidence interval percent to use for  
#'   lower (\code{cilow}) and upper (\code{ciup}) limits of bootrapped 
#'   quantiles. 
#'
#' @return list. Contains ...
#' @export
#'
#' @examples
#' load("data/ssb.Rdata")
#' source("R/qboot.R")
#' source("R/qbootplot.R")
#' X <- array(ssb, dim = dim(ssb)[c(2,6)])
#' dimnames(X) <- dimnames(ssb)[c(2,6)]
#' X <- t(X)
#' res <- qboot(X, nboot = 10, rerrTarg = 1)
#' qbootplot(res)
#' 
#' 
#' 
qboot <- function(X, nboot = 30, quant = c(0.95,0.75,0.5,0.25,0.05), 
  niter = exp(seq(log(nrow(X))*0.5, log(nrow(X)), length.out = 10)),
  estfun = "median", aggfun = "max", rerrTarg = 1, ci = 0.95, verbose = TRUE
){

  # add some tests of X dims
  X <- as.array(X)
  
  # empty matrix for bootstrapped quantiles
  qn <- array(NaN, dim=c(ncol(X), length(quant), nboot),
    dimnames = list(dimnames(X)[[2]], quant=quant, boot=seq(nboot)))
  # empty matrix for std.err and mean results
  est <- mad <- cilow <- ciup <- array(NaN, dim=c(ncol(X), length(quant), length(niter)),
    dimnames = list(dimnames(X)[[2]], quant=quant, niter=niter))
  for(i in seq(niter)){
    # determine the quantiles of sub-sampled iterations (bootstrapping)
    for(n in seq(nboot)){
      if(ncol(X)==1){
        qn[,,n] <- quantile(X[sample(nrow(X), niter[i], replace = TRUE),], prob=quant)
      } else {
        qn[,,n] <- t(apply(X[sample(nrow(X), niter[i], replace = TRUE),], 
          MARGIN = 2, FUN = quantile, prob=quant))
      }
    }
    # determine the CIs of bootstrapped quantiles
    Q <- apply(qn, MARGIN = 1:2, FUN = quantile, prob = c((1-ci)/2, 0.5, ci+(1-ci)/2))
    cilow[,,i] <- Q[1,,]
    est[,,i] <- Q[2,,]
    ciup[,,i] <- Q[3,,]
    # median absolute deviation of bootstrapped quantiles
    mad[,,i] <- apply(qn, MARGIN = 1:2, FUN = function(x){
      median(abs(x - median(x)))
    })
    if(verbose) print(paste0("niter = ", i, " of ", length(niter), " is finished"))
  }
  
  rerr <- mad / est * 100 # relative absolute error
  # rerr <- (ciup - cilow) / est * 100 # relative error 

  
  # fit log-log linear regression; lm(log(niter) ~ log(rerr) + quant)
  df <- expand.grid(dimnames(rerr))
  df$rerr <- c(rerr)
  agg <- aggregate(rerr ~ quant + niter, data = df, FUN = aggfun)
  agg$niter <- as.numeric(as.character(agg$niter))
  fit <- lm(log(niter) ~ log(rerr) + quant, agg)
  newdat <- data.frame(rerr = rerrTarg, quant = factor(quant))
  newdat$niter <- exp(predict(fit, newdata = newdat))
  
  # return results
  ret <- list(est = est, cilow = cilow, ciup = ciup, mad = mad, rerr = rerr,
    fit = fit, nsim = newdat,
    nboot = nboot, quant = quant, niter = niter, agg = agg,
    aggfun = aggfun, estfun = estfun, ci = ci, rerrTarg = rerrTarg)
  class(ret) <- "qboot"
  return(ret)
}
