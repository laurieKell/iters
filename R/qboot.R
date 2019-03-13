#' Calculate MSE quantile estimation error via bootstrapping and predict
#' number of iterations needed for a given precision
#'
#' @param X matrix. Statistic of interest across years (rows) and 
#'   iterations (columns)
#' @param nboot numeric. Number of bootstraps. 
#' @param quant vector. Quantile probabilities 
#'   (Default: \code{quant = c(0.95,0.75,0.5,0.25,0.05)}).
#' @param niter vector. Cumulative number of iterations to perform bootstrapped
#'   estimates (Default: \code{niter = exp(seq(log(ncol(X))*0.5, 
#'   log(ncol(X)), length.out = 20))}).
#' @param aggfun string. Name of function to apply to quantile relative error by year
#'   for use in log-log model (Default: \code{aggfun = "max"}). Use of "max" is 
#'   consistent with 'Prob3', which is the maximum annual probability of SSB 
#'   dropping below Blim. 'Prob2' would be the average annual probability 
#'   of SSB dropping below Blim (i.e. \code{aggfun = "mean"}).    
#' @param targsemu numeric. Percent relative error (se/mu) desired.
#' @param verbose logical. Should progress be printed 
#'   (Default: code\{verbose = TRUE})
#'
#' @return list. Contains 
#' @export
#'
#' @examples
#' load("data/SSB.Rdata")
#' X <- array(ssb, dim = dim(ssb)[c(2,6)])
#' dimnames(X) <- dimnames(ssb)[c(2,6)]
#' res <- qboot(X, nboot = 10)
#' 
#' 
#' 
qboot <- function(X, nboot = 20, quant = c(0.95,0.75,0.5,0.25,0.05), 
  niter = exp(seq(log(ncol(X))*0.5, log(ncol(X)), length.out = 20)),
  aggfun = "max", targsemu = 1, verbose = TRUE
){

  # add some tests of X dims
  
  # empty matrix for bootstrapped quantiles
  Qn <- array(NaN, dim=c(nrow(X), length(quant), nboot),
    dimnames = list(year=seq(nrow(X)), quant=quant, boot=seq(nboot)))
  # empty matrix for std.err and mean results
  SE <- MU <- array(NaN, dim=c(nrow(X), length(quant), length(niter)),
    dimnames = list(year=seq(nrow(X)), quant=quant, niter=niter))
  for(i in seq(niter)){
    for(n in seq(nboot)){
      if(nrow(X)==1){
        Qn[,,n] <- quantile(X[,sample(ncol(X), niter[i], replace = TRUE)], prob=quant)
      } else {
        Qn[,,n] <- t(apply(X[,sample(ncol(X), niter[i], replace = TRUE)], 
          MARGIN = 1, FUN = quantile, prob=quant))
      }
    }
    MU[,,i] <- apply(Qn, MARGIN = 1:2, FUN = mean)
    SE[,,i] <- apply(Qn, MARGIN = 1:2, FUN = sd)

    if(verbose) print(paste0("niter = ", i, " of ", length(niter), " is finished"))
  }
  
  SEMU <- SE/MU*100
  
  # fit log-log linear regression; lm(log(niter) ~ log(semu) + quant)
  df <- expand.grid(dimnames(SEMU))
  df$semu <- c(SEMU)
  agg <- aggregate(semu ~ quant + niter, data = df, FUN = aggfun)
  agg$niter <- as.numeric(as.character(agg$niter))
  fit <- lm(log(niter) ~ log(semu) + quant, agg)
  newdat <- data.frame(semu = targsemu, quant = factor(quant))
  newdat$niter <- exp(predict(fit, newdata = newdat))
  
  # return results
  ret <- list(SE = SE, MU = MU, SEMU = SEMU, fit = fit, nsim = newdat,
    nboot = nboot, quant = quant, niter = niter, agg = agg,
    aggfun = aggfun, targsemu = targsemu)
  class(ret) <- "qboot"
  return(ret)
}
