#############################################
####### Simulation: Uniform #################
#############################################
# Author: Matias D. Cattaneo, Max Farrell, Yingjie Feng
# Date:   05/19/2019

rm(list=ls())
setwd("E:/Dropbox/Bspline/Simulations/replication")
# load package
library(lspartition)
library(matrixStats)
library(foreach)
library(doParallel)
library(doRNG)
library(Hmisc)
source("preamble_rep.R")

########## Uniform Inference ##########
# pars
smooth   <- NULL
par <- read.csv("model.csv", header = T, colClasses=c(rep("numeric", 5), rep("character", 2),
                                                      "logical", "numeric", "numeric"))

##par value###################
num <- c(1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14, 15:19)
for (ell in 1:19) {
  j <- num[ell]
  optk  <- as.matrix(read.table(paste("post/raw/optkappa_", j, ".txt", sep=""), header = F))
  range <- min(optk):max(optk)
  model <- funlist[[par$model[j]]]
  sigma <- par$sigma[j]
  n     <- par$n[j]
  d     <- par$d[j]
  m     <- par$m[j]
  q     <- par$q[j]
  method  <- par$method[j]
  ktype   <- par$ktype[j]
  proj    <- par$proj[j]
  kimse   <- round(par$kopt[j])
  rep <- 5000

  if (d == 1) {
    grid <- as.matrix(seq(0, 1, length.out = 102)[-c(1, 102)])
    tval <- model(grid[,1])
  }

  if (d == 2) {
    grid <- seq(0, 1, length.out = 12)[-c(1,12)]
    grid <- expand.grid(grid, grid)
    tval <- model(grid[,1], grid[,2])
  }

  if (d == 3) {
    grid  <- seq(0, 1, length.out = 7)[-c(1,7)] + 0.05 # just to avoid always evaluating at knots
    gridz <- seq(0, 1, length.out = 6)[-c(1,6)]
    grid  <- expand.grid(grid, grid, gridz)
    tval  <- model(grid[,1], grid[,2], grid[,3])
  }

  ## Uniform band simulation ########
  cl <- makeCluster(3)
  registerDoParallel(cl)
  writeLines(c(""), "log.txt")

  band <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('lspartition'),
                  .combine=rbind) %dorng% {
                   sink("log.txt", append=TRUE)
                   cat(paste("Starting iteration", i, "\n"))
                   sink()
                   output <- sim.unicov2(i, n, d, model, grid, method, m, q, smooth, ktype,
                                         kimse, proj, tval, optk, sigma)
                   output
                 }

  stopCluster(cl)


  # Extract information###############
  write.table(band, paste("post/raw/final.rej.full", j, "txt", sep = "."), row.names = F, col.names = F)

  #band <- read.table(paste("post/raw/final.rej.full",j,"txt", sep="."))
  if ((j<=14 & ell%%2!=0) | j>14) table <- matrix(NA, nrow = 24, ncol = 8)
  ncol <- 3; L <- nrow(grid); gap <- 8 * L + 8
  rowname <- rep(1:gap, rep)

  rejlist <- list()

  for (i in 1:8) {
    rejlist[[i]] <- cbind(rep(1:L, rep), band[rowname %in% ((i-1)*L+1):(i*L),])
  }

  length <- list()
  for (i in 1:8) {
    length[[i]] <- colMeans(band[rowname == 8*L+i,])
  }
  names(length) <- c("pl.cl",  "wb.cl",
                     "pl.bc1", "wb.bc1",
                     "pl.bc2", "wb.bc2",
                     "pl.bc3", "wb.bc3")

  bandcov <- lapply(rejlist, function(x) banderror(x, L, rep))
  names(bandcov) <- c("pl.cl",  "wb.cl",
                      "pl.bc1", "wb.bc1",
                      "pl.bc2", "wb.bc2",
                      "pl.bc3", "wb.bc3")

  ##### Horowitz measures ######################
  for (i in 1:8) {
    rejlist[[i]] <- rejlist[[i]][order(rejlist[[i]][,1]),]
  }
  griderrors <- lapply(rejlist, function(x) pcoverror(x, L, rep))
  proportion <- lapply(griderrors, function(x) cpr(x[,-1]))
  names(proportion) <- c("pl.cl",  "wb.cl",
                         "pl.bc1", "wb.bc1",
                         "pl.bc2", "wb.bc2",
                         "pl.bc3", "wb.bc3")

  aveerror   <- lapply(griderrors, function(x) colMeans(x[,-1]))
  names(aveerror) <- c("pl.cl",  "wb.cl",
                       "pl.bc1", "wb.bc1",
                       "pl.bc2", "wb.bc2",
                       "pl.bc3", "wb.bc3")

  if (ktype == "uni") {
    table[,1:4] <- rbind(cbind(proportion$pl.cl,  aveerror$pl.cl,  length$pl.cl,  bandcov$pl.cl),
                         cbind(proportion$pl.bc1, aveerror$pl.bc1, length$pl.bc1, bandcov$pl.bc1),
                         cbind(proportion$pl.bc2, aveerror$pl.bc2, length$pl.bc2, bandcov$pl.bc2),
                         cbind(proportion$pl.bc3, aveerror$pl.bc3, length$pl.bc3, bandcov$pl.bc3),
                         cbind(proportion$wb.cl,  aveerror$wb.cl,  length$wb.cl,  bandcov$wb.cl),
                         cbind(proportion$wb.bc1, aveerror$wb.bc1, length$wb.bc1, bandcov$wb.bc1),
                         cbind(proportion$wb.bc2, aveerror$wb.bc2, length$wb.bc2, bandcov$wb.bc2),
                         cbind(proportion$wb.bc3, aveerror$wb.bc3, length$wb.bc3, bandcov$wb.bc3))
  }

  if (ktype == "qua") {
    table[,5:8] <- rbind(cbind(proportion$pl.cl,  aveerror$pl.cl,  length$pl.cl,  bandcov$pl.cl),
                         cbind(proportion$pl.bc1, aveerror$pl.bc1, length$pl.bc1, bandcov$pl.bc1),
                         cbind(proportion$pl.bc2, aveerror$pl.bc2, length$pl.bc2, bandcov$pl.bc2),
                         cbind(proportion$pl.bc3, aveerror$pl.bc3, length$pl.bc3, bandcov$pl.bc3),
                         cbind(proportion$wb.cl,  aveerror$wb.cl,  length$wb.cl,  bandcov$wb.cl),
                         cbind(proportion$wb.bc1, aveerror$wb.bc1, length$wb.bc1, bandcov$wb.bc1),
                         cbind(proportion$wb.bc2, aveerror$wb.bc2, length$wb.bc2, bandcov$wb.bc2),
                         cbind(proportion$wb.bc3, aveerror$wb.bc3, length$wb.bc3, bandcov$wb.bc3))
  }
  table[,c(2,3,6,7)] <- round(table[,c(2,3,6,7)], 3)
  table[,c(4,8)]     <- round(table[,c(4,8)], 1)
  
}






#######################################################
# Plot CB
j <- 1
optk  <- as.matrix(read.table(paste("output/raw/optkappa_", j, ".txt", sep=""), header = F))
range <- min(optk):max(optk)
model <- funlist[[par$model[j]]]
sigma    <- par$sigma[j]
n <- par$n[j]
d <- par$d[j]
m        <- par$m[j]
q        <- par$q[j]
method   <- par$method[j]
ktype    <- par$ktype[j]
proj     <- par$proj[j]
kimse    <- round(par$kopt[j])
if (d == 1) {
  grid <- as.matrix(seq(0, 1, length.out = 102)[-c(1, 102)])
}
L <- nrow(grid)

# Comparison: one realization
set.seed(1234)
data  <- dgp(n, model, d, sigma)
dev.new(width=5, height=3)
pdf(file="E:/Dropbox/Bspline/Simulations/replication/post/fig_cb.pdf") 
plot(grid[,1], model(grid[,1]), type="l", lty=1, col=1, 
     xlab=c("x"), ylab=c("y"),ylim=c(-2,1),cex.axis=1.5, cex.lab=1.6)
k.dpi <- lspkselect(data$y, data$x, m=m, method=method, vce="hc2",
                    ktype=ktype, kselect="all", rotnorm=F)$ks[,c("k.dpi")]
band <- plotband(data, d, grid, method, m, q, J, smooth, ktype, k.dpi, proj)

# point estimate
lines(grid[,1], band[1:L], lty=2, col=2)
ll <- seq(2,16,2)
for (i in 1:length(ll)) {
  j <- ll[i]
  lines(grid[,1], band[((j-1)*L+1):(j*L)], lty=i+2, col=i+2)
  lines(grid[,1], band[((j*L)+1):((j+1)*L)], lty=i+2, col=i+2)
}
legend("bottomright",legend=c("True", "Est.", "P.0", "W.0",
                              "P.1", "W.1","P.2", "W.2","P.3", "W.3"),
       lty=1:10, col=1:10, ncol=5)
dev.off()

