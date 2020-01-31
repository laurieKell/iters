
load(file = "data/codMSE_A.Rdata")
source("R/pcalc.R")
source("R/qboot.R")
source("R/qbootplot.R")

M <- t(A[1,1,1,,]*0)
dim(M)

df <- expand.grid(dimnames(A)[5:4], stringsAsFactors = FALSE)
df$niter <- 0
df$risk <- NaN
df$medC <- NaN
df$medSSB <- NaN
df$varC <- NaN
df$cont <- TRUE
df$medC.rank <- 0
df$q05ciupp <- NaN
df$rerr <- NaN

iter.min <- 100
iter.max <- 1000
iter.step <- 100
rerrTarg <- 2 # target relative error (as %)
years <- tail(dimnames(A)$year, 10)
Blim <- 107000
set.seed(1)

i <- iter.min
while(i <= iter.max & sum(df$cont) > 1){
  # i = 1; j = 1
  # run more iterations
  df$niter[which(df$cont)] <- df$niter[which(df$cont)] + iter.step
  print(paste(sum(df$cont), "grids continuing to niter", i ))
  
  pb <- txtProgressBar(min = 1, max = nrow(df), style = 3)
  for(j in seq(nrow(df))){
    if(df$cont[j]){
      X <- A["SSB", years, as.character(seq(df$niter[j])), df$Ftrgt[j], df$Btrig[j]]
      df$risk[j] <- pcalc(SSB = X, Blim = Blim)$P3
      df$medC[j] <- median(A["Catch", years, as.character(seq(df$niter[j])), df$Ftrgt[j], df$Btrig[j]])
      df$medSSB[j] <- median(A["SSB", years, as.character(seq(df$niter[j])), df$Ftrgt[j], df$Btrig[j]])
      df$varC[j] <- median(apply(A["Catch", years, as.character(seq(df$niter[j])), df$Ftrgt[j], df$Btrig[j]], 2, function(x){sd(log(x))}))
      # df[j,]
      # matplot(X, col = adjustcolor(1,0.1), t = "l", lty = 1, pch = NA, ylim = c(0,max(X)) ); abline(h = Blim, lty = 2, col = 2)
      qres <- qboot(X = t(X), nboot = 20, quant = c(0.05, 0.5, 0.95), 
        niter = exp(seq(log(ncol(X))*0.5, log(ncol(X)), length.out = 3)),
        aggfun = "max", rerrTarg = rerrTarg, ci = 0.95, verbose = FALSE)
      df$q05ciupp[j] <- min(qres$ciup[, "0.05", as.character(df$niter[j])])
      # df$q05cilow[j] <- min(qres$cilow[, "0.05", as.character(df$niter[j])])
      
      df$rerr[j] <- max(qres$rerr[, "0.05", as.character(df$niter[j])])
      # max(qres$rerr[, "0.05", 1])
      # qbootplot(qres)
    }
    # print(paste("niter =", i, "| grid =", j))
    setTxtProgressBar(pb, j)
  }
  close(pb)

  ### criteria for continuing with simulations
  
  # stop risky ones
  # df$cont <- TRUE
  rmv <- which(df$q05ciupp < Blim)
  if(length(rmv)>0) df$cont[rmv] <- FALSE

  # stop if not in top X% 
  df$medC.rank <- NaN
  df$medC.rank[which(df$cont)] <- rank(df$medC[which(df$cont)])
  maxRank <- max(df$medC.rank, na.rm = TRUE)
  df$medC.rank <- maxRank - df$medC.rank + 1
  
  # keep top 75% or in top 30
  rankThresh <- max(c(0.75*maxRank), 30)
  rmv <- which(df$medC.rank > rankThresh)
  if(length(rmv)>0) df$cont[rmv] <- FALSE

  # stop when rerrTarg is reached for all continuing
  if(all(df$rerr[df$cont] < rerrTarg)) df$cont <- FALSE
  
  m <- M; m[] <- df$cont; image(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt), z = m)
 
  i <- i + iter.step
}



# Best
fitness <- df$medC * (df$risk < 0.05) * (df$niter==max(df$niter))
ra <- rank(fitness)
top <- df[which((max(ra) - ra + 1) %in% 1:10),]
top <- top[order(top$medC, decreasing = T),]
top
# should be Btrig = 170000, Ftrgt = 0.38
actualBest <- df[which(df$Btrig=="170000" & df$Ftrgt == "0.38"),]
actualBest

# fraction of max iterations needed
sum(df$niter) / (nrow(df)*iter.max)


B <- t(A[1,1,1,,]) * NaN
x = as.numeric(dimnames(B)$Btrig)
y = as.numeric(dimnames(B)$Ftrgt)

RISK <- B
RISK[] <- df$risk
image(x = x, y = y, z = RISK)
contour(x = x, y = y, z = RISK, add = TRUE)
points(as.numeric(top$Btrig[1]), as.numeric(top$Ftrgt[1]), col=4)
points(as.numeric(actualBest$Btrig[1]), as.numeric(actualBest$Ftrgt[1]), col=3, pch=20)

MEDC <- B
MEDC[] <- df$medC
image(x = x, y = y, z = MEDC, col = pals::parula(100))
contour(x = x, y = y, z = MEDC, add = TRUE)
contour(x = x, y = y, z = RISK, add = TRUE, col = 2)#, levels = 0.05)
points(as.numeric(top$Btrig[1]), as.numeric(top$Ftrgt[1]), col=4)
points(as.numeric(actualBest$Btrig[1]), as.numeric(actualBest$Ftrgt[1]), col=3, pch=20)

MEDSSB <- B
MEDSSB[] <- df$medSSB
image(x = x, y = y, z = MEDSSB, col = pals::parula(100))
contour(x = x, y = y, z = MEDSSB, add = TRUE)
points(as.numeric(top$Btrig[1]), as.numeric(top$Ftrgt[1]), col=4)
points(as.numeric(actualBest$Btrig[1]), as.numeric(actualBest$Ftrgt[1]), col=3, pch=20)

NITER <- B
NITER[] <- df$niter
image(x = x, y = y, z = NITER)
contour(x = x, y = y, z = NITER, add = TRUE)
contour(x = x, y = y, z = RISK, add = TRUE, col = 2)#, levels = 0.05)
points(as.numeric(top$Btrig[1]), as.numeric(top$Ftrgt[1]), col=4)
points(as.numeric(actualBest$Btrig[1]), as.numeric(actualBest$Ftrgt[1]), col=3, pch=20)



