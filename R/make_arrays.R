# devtools::install_github("flr/mse")

library(FLCore)
library(mse)
library(pals)

datPath <- "data/codMSE/extracted"
L <- list.files(datPath)
stk <- readRDS(file = file.path(datPath, L[[1]]) )

L <- L[-which(L == "stats.rds")]
L
length(L)



S <- strsplit(x = L, fixed = T, split = "_")
fun <- function(x) {
  data.frame(stk = x[1], hcr = x[2], 
    Ftrgt = strsplit(x = x[3], split = "-")[[1]][2],
    Btrig = strsplit(x = x[4], split = "-")[[1]][2],
    TACconstr = strsplit(x = x[5], split = "-")[[1]][2], 
    BB = strsplit(strsplit(x = x[6], split = "-")[[1]][2], ".", fixed = T)[[1]][1], 
    stringsAsFactors = FALSE)
}
D <- do.call("rbind", lapply(S, FUN = fun))
D$fname = L
str(D)

A <- array(NaN, 
  dim = c(2, dim(stk@stock)[2], dim(stk@stock)[6], length(unique(D$Ftrgt)), 
    length(unique(D$Btrig))), 
  dimnames = list(
    data = c("SSB", "Catch"), 
    year = dimnames(stk@stock)[2][[1]],
    iter = dimnames(stk@stock)[6][[1]],
    Ftrgt = as.character(sort(as.numeric(unique(D$Ftrgt)))),
    Btrig = as.character(sort(as.numeric(unique(D$Btrig))))
  )
)
  
dim(A)
dimnames(A)

for(i in seq(L)){
  stk.i <- readRDS(file = file.path(datPath, L[[i]]))
  A["SSB",,,as.character(D$Ftrgt[i]),as.character(D$Btrig[i])] <- array(ssb(stk.i@stock), dim = dim(stk.i@stock)[c(2,6)])
  A["Catch",,,as.character(D$Ftrgt[i]),as.character(D$Btrig[i])] <- array(catch(stk.i@stock), dim = dim(stk.i@stock)[c(2,6)])
  print(paste(i, "is finished"))
}
beepr::beep()


image(A[1,,,1,1], useRaster = T)
image(A[1,,,10,8], useRaster = T)

matplot(x = as.numeric(dimnames(A)$year), y = A[1,,,1,1], col = adjustcolor(1,0.1), t = "l", lty = 1, pch = NA)
matplot(x = as.numeric(dimnames(A)$year), y = A[1,,,40,11], col = adjustcolor(1,0.1), t = "l", lty = 1, pch = NA)


save(A, file = "data/codMSE_A.Rdata")
# load(file = "data/codMSE_A.Rdata")




# risk of SSB < Blim over last 10 years of sim

prisk <- expand.grid(Btrig = dimnames(A)$Btrig, Ftrgt = dimnames(A)$Ftrgt, 
  P1 = NaN, P2 = NaN, P3 = NaN, stringsAsFactors = FALSE)
prisk
dim(prisk)
head(prisk)

source("R/pcalc.R")

Ftrgt <- 0.5
Btrgr <- 140000
Blim <- 107000
SB <- A["SSB",,, as.character(Ftrgt), as.character(Btrgr)]
pcalc(SSB = SB[tail(seq(dim(SB)[1]), 10),], Blim = Blim)
matplot(x = as.numeric(dimnames(A)$year), y = SB, col = adjustcolor(1,0.1), t = "l", lty = 1, pch = NA)
abline(h = Blim, lty = 2, col = 2)

pb <- txtProgressBar(min = 1, max = nrow(prisk), style = 3)
for(i in seq(nrow(prisk))){
  SB <- A["SSB",,, as.character(prisk$Ftrgt[i]), as.character(prisk$Btrig[i])]
  prisk[i, c("P1", "P2", "P3")] <- pcalc(SSB = SB[tail(seq(dim(SB)[1]), 10),], 
    Blim = Blim )
  setTxtProgressBar(pb, i)
}
close(pb)

P <- A[1,1,1,,] * NaN
dim(P)

P <- array(prisk$P3, dim = dim(A)[5:4])
image(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt),
  z = P, col = parula(100))
contour(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt),
  z = P, add = T, col = "white")



dimnames(A)

medC <- t(apply(A["Catch", tail(seq(dim(A)[2]), 10),,,], MARGIN = 3:4, FUN = median))
image(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt), 
  z = medC, col = rev(parula(50)))
contour(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt),
  z = medC, add = T, col = "white")


sumC <- t(apply(A["Catch", tail(seq(dim(A)[2]), 10),,,], MARGIN = 3:4, FUN = sum))
image(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt), 
  z = sumC, col = rev(parula(50)))
contour(x = as.numeric(dimnames(A)$Btrig), y = as.numeric(dimnames(A)$Ftrgt),
  z = sumC, add = T, col = "white")

