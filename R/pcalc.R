pcalc <- function(SSB, Blim){
  # annual probability of SSB < Blim
  annP <- apply(SSB, 1, function(x){mean(x < Blim)})
  # Ps
  P1 <- mean(annP)
  P2 <- mean(apply(SSB, 2, function(x){max(x < Blim)}))
  P3 <- max(annP)
  return(list(P1 = P1, P2 = P2, P3 = P3))
}
