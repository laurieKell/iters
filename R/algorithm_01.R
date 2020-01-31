
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
rerrTarg <- 1 # target relative error (as %)
years <- tail(dimnames(A)$year, 10)
Blim <- 107000
set.seed(1)

doplot <- TRUE
B <- t(A[1,1,1,,]) * NaN
x = as.numeric(dimnames(B)$Btrig)
y = as.numeric(dimnames(B)$Ftrgt)


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
      # matplot(X, col = adjustcolor(1,0.1), t = "l", lty = 1, pch = NA, ylim = c(0,max(X)) ); abline(h = Blim, lty = 2, col = 2)
      qres <- qboot(X = t(X), nboot = 20, quant = c(0.05, 0.5, 0.95), 
        niter = ncol(X), ci = 0.95, verbose = FALSE)
      df$q05ciupp[j] <- min(qres$ciup[, "0.05", as.character(df$niter[j])])
      df$rerr[j] <- max(qres$rerr[, "0.05", as.character(df$niter[j])])

    }
    setTxtProgressBar(pb, j)
  }
  close(pb)

  ### criteria for continuing with simulations
  
  # stop risky ones (if upper CI of 5% quant is below Blim)
  rmv <- which(df$q05ciupp < Blim)
  if(length(rmv)>0) df$cont[rmv] <- FALSE

  # stop if not in top X% 
  df$medC.rank <- NaN
  df$medC.rank[which(df$cont)] <- rank(df$medC[which(df$cont)])
  maxRank <- max(df$medC.rank, na.rm = TRUE)
  df$medC.rank <- maxRank - df$medC.rank + 1
  
  # keep top 75% or minimum of top 30
  rankThresh <- max(c(0.75*maxRank), 30)
  rmv <- which(df$medC.rank > rankThresh)
  if(length(rmv)>0) df$cont[rmv] <- FALSE

  # stop when rerrTarg is reached for all continuing
  if(all(df$rerr[df$cont] < rerrTarg)) df$cont <- FALSE
  
  # plot
  if(doplot){
    png(file.path("figs", paste0("algo_iter_", sprintf(fmt = "%04d", i), ".png")), width = 5, height = 5, 
      res = 400, units = "in")
    op <- par(mfrow = c(2,2), mar = c(3.5,3.5,1.5,0.5), mgp = c(2,0.5,0), ps = 10, oma = c(0,0,1.5,0))
    # risk
    risk <- B; risk[] <- df$risk
    
    # best
    best <- df[which.max(df$cont * (df$risk < 0.05) * df$medC),]

    # medC
    medC <- B; medC[] <- df$medC
    image(x = x, y = y, z = medC, col = grey.colors(100),
      xlab = "Btrigger", ylab = "Ftarget")
    points(as.numeric(best$Btrig), as.numeric(best$Ftrgt), col = 3, pch = 16)
    contour(x = x, y = y, z = medC, col = 1, add = TRUE)
    contour(x = x, y = y, z = risk, col = 2, add = TRUE)
    mtext(text = "median Catch", side = 3, line = 0.25, adj = 0)
    
    # medSSB
    medSSB <- B; medSSB[] <- df$medSSB
    image(x = x, y = y, z = medSSB, col = grey.colors(100),
      xlab = "Btrigger", ylab = "Ftarget")
    points(as.numeric(best$Btrig), as.numeric(best$Ftrgt), col = 3, pch = 16)
    contour(x = x, y = y, z = medSSB, col = 1, add = TRUE)
    contour(x = x, y = y, z = risk, col = 2, add = TRUE)
    mtext(text = "median SSB", side = 3, line = 0.25, adj = 0)
  
    # iterations
    niter <- B; niter[] <- df$niter
    image(x = x, y = y, z = niter, col = grey.colors(100),
      xlab = "Btrigger", ylab = "Ftarget")
    points(as.numeric(best$Btrig), as.numeric(best$Ftrgt), col = 3, pch = 16)
    contour(x = x, y = y, z = niter, col = 1, add = TRUE)
    contour(x = x, y = y, z = risk, col = 2, add = TRUE)
    mtext(text = "Iterations completed", side = 3, line = 0.25, adj = 0)
    
    # cont
    cont <- B; cont[] <- df$cont
    image(x = x, y = y, z = cont, col = grey.colors(100),
      xlab = "Btrigger", ylab = "Ftarget")
    points(as.numeric(best$Btrig), as.numeric(best$Ftrgt), col = 3, pch = 16)
    contour(x = x, y = y, z = risk, col = 2, add = TRUE)
    mtext(text = "Continuing", side = 3, line = 0.25, adj = 0)
    
    mtext(text = paste("iter =", i), side = 3, line = 0.25, adj = 0.1, outer = TRUE)
    par(op)
    dev.off()
  }
  
  i <- i + iter.step
}

best


library(magick)
fnames <- list.files(path = "figs", pattern = "iter")
for(i in seq(fnames)){
  if(i == 1)  img <- image_read(file.path("figs", fnames[i]))
  if(i > 1){
    frame.i <- image_read(file.path("figs", fnames[i]))
    img <- append(img, frame.i)
  }
}

img <- image_scale(img, geometry = geometry_size_percent(width = 50, height = NULL))
img <- image_animate(img, optimize = TRUE, delay = 100)
image_info(img)

image_write(img, path = "figs/algo_animation.gif", format = "gif")





bigdata <- image_read('https://jeroen.github.io/images/bigdata.jpg')
frink <- image_read("https://jeroen.github.io/images/frink.png")
logo <- image_read("https://jeroen.github.io/images/Rlogo.png")
img <- c(bigdata, logo, frink)
img <- image_scale(frame.i, "300x300")
image_info(img)

animation <- image_animate(img, fps = 2, optimize = TRUE)
print(animation)

system('"C:\\Program Files\\ImageMagick-7.0.7-Q16\\convert.exe" -delay 20 -loop 0 figs/algo_iter*.png figs/algo_iter_animation.gif')
system('"C:\\Program Files\\ImageMagick-7.0.7-Q16\\convert.exe" -gravity Center -crop 2800x1200+0+0 -delay 20 -loop 0 figs/algo_iter*.png figs/algo_iter_animation.gif')


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
contour(x = x, y = y, z = RISK, add = TRUE, col = 2)
points(as.numeric(top$Btrig[1]), as.numeric(top$Ftrgt[1]), col=4)
points(as.numeric(actualBest$Btrig[1]), as.numeric(actualBest$Ftrgt[1]), col=3, pch=20)



