#' Calculate MSE quantile estimation error via bootstrapping and predict
#' number of iterations needed for a given precision.
#' 
#' @description Calculates bootstrapped estimates of quantiles (median). 
#'   Estimation error (median absolute deviation (MAD)) and confidence
#'   intervals (CI) are also estimated. The summary statistic of 
#'   relative error (rerr=(MAD/est)*100) is also returned. When estimates are 
#'   performed over various levels of iterations (`niter`), a linear 
#'   regression is fit to predict the number of iterations needed to achieve
#'   a given level of error (`rerrTarget`).
#' 
#' @param X matrix. Statistic of interest across years (columns) and 
#'   iterations (rows)
#' @param nboot numeric. Number of bootstraps. 
#' @param quant vector. Quantile probabilities 
#'   (Default: \code{quant = c(0.95,0.75,0.5,0.25,0.05)}).
#' @param niter vector. Cumulative number of iterations to perform bootstrapped
#'   estimates (Default: \code{niter = exp(seq(log(nrow(X))*0.5, log(nrow(X)), 
#'   length.out = 5))}).
#' @param ci numeric. Confidence interval level to use for  
#'   lower (\code{cilow}) and upper (\code{ciup}) limits of bootrapped 
#'   quantiles (Default: `ci = 0.95`). 
#' @param aggfun string. Name of function to apply to quantile relative error 
#'   across years (Default: \code{aggfun = "max"}). Use of "max" is 
#'   consistent with 'Prob3', which is the maximum annual probability of SSB 
#'   dropping below Blim. 'Prob2' would be the average annual probability 
#'   of SSB dropping below Blim (i.e. \code{aggfun = "mean"}). Statistic is 
#'   applied to each level of `niter``.`   
#' @param rerrTarg numeric. Desired relative error (as a percent; Default: 1)
#'   for which to predict the required number of iterations. Only applicable if
#'   three or more levels of `niter` are bootstrapped, allowing for the fitting 
#'   (and prediction) of relative error as a function of number of iterations.
#' @param verbose logical. Should progress be printed 
#'   (Default: code\{verbose = TRUE})
#'   
#'   
#' @return list. Contains bootstrapped estimates of quantiles (`quant`) 
#'   (estimate: `est`)
#' @export
#'
#' @examples
#' load("data/ssb.Rdata")
#' source("R/qboot.R")
#' source("R/qbootplot.R")
#' X <- array(ssb, dim = dim(ssb)[c(2,6)])
#' dimnames(X) <- dimnames(ssb)[c(2,6)]
#' X <- t(X)
#' res <- qboot(X, nboot = 100, rerrTarg = 1)
#' res$ciLevs
#' qbootplot(res)
#' 
#' res <- qboot(X, nboot = 100, rerrTarg = 1, niter = 200)
#' # qbootplot(res) # will return an error as there is nothing to plot
#' 
#' 
#' 
qboot <- function(X, nboot = 30, quant = c(0.95,0.75,0.5,0.25,0.05), 
  niter = exp(seq(log(nrow(X))*0.5, log(nrow(X)), length.out = 5)),
  ci = 0.95, aggfun = "max", rerrTarg = 1, verbose = TRUE
){

  # add some tests of X dims
  X <- as.array(X)
  
  # empty matrix for bootstrapped quantiles
  qn <- array(NaN, dim=c(ncol(X), length(quant), nboot),
    dimnames = list(dimnames(X)[[2]], quant=quant, boot=seq(nboot)))
  # empty matrix for std.err and mean results
  est <- mad <- cilow <- ciup <- array(NaN, dim=c(ncol(X), length(quant), length(niter)),
    dimnames = list(dimnames(X)[[2]], quant=quant, niter=niter))
  
  if(verbose & length(niter) > 1) pb <- txtProgressBar(min = 1, max = length(niter), style = 3)
  for(i in seq(length(niter))){
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
    if(verbose & length(niter) > 1) setTxtProgressBar(pb, i)
  }
  if(verbose & length(niter) > 1) close(pb)

  # relative error of MAD
  rerr <- mad / est * 100 
  
  # fit log-log linear regression; lm(log(niter) ~ log(rerr) + quant)
  if(length(niter) >= 3){
    df <- expand.grid(dimnames(rerr))
    df$rerr <- c(rerr)
    agg <- aggregate(rerr ~ quant + niter, data = df, FUN = aggfun)
    agg$niter <- as.numeric(as.character(agg$niter))
    fit <- lm(log(niter) ~ log(rerr) + quant, agg)
    newdat <- data.frame(rerr = rerrTarg, quant = factor(quant))
    newdat$niter <- exp(predict(fit, newdata = newdat))
  }else{
    if(verbose){
      print("length(niter) < 3; No fitting of model conducted: log(niter)~log(rerr)+quant")
    }
    agg <- NULL
    fit = NULL
    newdat = NULL
  }
  
  # return results
  ret <- list(
    est = est, cilow = cilow, ciup = ciup, mad = mad, rerr = rerr,
    fit = fit, fitpred = newdat,
    nboot = nboot, quant = quant, niter = niter, agg = agg,
    aggfun = aggfun, ci = ci, ciLevs = c((1-ci)/2, ci+(1-ci)/2), rerrTarg = rerrTarg)
  class(ret) <- "qboot"
  return(ret)
}
