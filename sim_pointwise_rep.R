#############################################
####### Simulation: Pointwise ###############
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

########## Pointwise Inference #######

### Pars ####
smooth   <- NULL
par <- read.csv("model.csv", header = T, colClasses=c(rep("numeric", 5), rep("character", 2),
                                                      "logical", "numeric", "numeric"))
for (j in 1:19) {

  sigma    <- par$sigma[j]
  kimse <- round(par$kopt[j])
  model <- funlist[[par$model[j]]]
  n <- par$n[j]
  d <- par$d[j]
  m        <- par$m[j]
  q        <- par$q[j]
  method   <- par$method[j]
  ktype    <- par$ktype[j]
  proj     <- par$proj[j]

  # range of kappa
  if (method=="bs") {
    range <- (kimse-2):(kimse+2)
    if (d == 3) range <- 1:3
  } else {
    if (d == 1) range <- 3:6
    if (d == 2) range <- 3:4
  }
  
  # repetitions
  rep <- 5000
  
  # three evaluating points
  if (d == 1) eval <- as.matrix(c(0.2, 0.5, 0.8))
  if (d == 2) eval <- rbind(c(0.5, 0.5), c(0.1, 0.5), c(0.1, 0.1))  
  if (d == 3) eval <- rbind(c(0.5, 0.5, 0.5), c(0.1, 0.1, 0.5), c(0.1, 0.1, 0.1))

  #########simulation: pointwise######
  ptm<-proc.time()
  cl <- makeCluster(3)
  registerDoParallel(cl)
  writeLines(c(""), "log.txt")

  output <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('lspartition'),
                     .combine=rbind) %dorng% {
                     sink("log.txt", append=TRUE)
                     cat(paste("Starting iteration", i, "\n"))
                     sink()
                     output <- sim(i, n, d, eval, method, m, q, smooth, ktype, range, 
                                   proj, model, sigma)
                     output
                    }

  stopCluster(cl)
  proc.time()-ptm
  
  write.table(output, paste("post/raw/final_pwtable", j, "txt", sep = "."), sep = ",", row.names = F, col.names = F)

  ####################
  # Save information #
  ####################
  if (d == 1) tval <- model(eval[,1])
  if (d == 2) tval <- model(eval[,1], eval[,2])
  if (d == 3) tval <- model(eval[,1], eval[,2], eval[,3])
  
  ncol <- length(range)+2; neval <- nrow(eval); gap <- neval * 16 + 1
  rowname <- rep(1:gap, rep)
  L <- neval * 4
  
  # Optimal kappa
  kappa <- output[1, 1:length(range)]
  rot   <- output[rowname == 1, length(range)+1]
  dpi   <- output[rowname == 1, length(range)+2]
  write.table(cbind(rot, dpi), paste("post/raw/optkappa_", j, ".txt", sep = ""), row.names = F, col.names = F)
  
  # Estimation and inference
  ptest <- output[rowname %in% (1+1):(1+L),] 
  se    <- output[rowname %in% (1+L+1):(1+2*L),]
  rej   <- output[rowname %in% (1+2*L+1):(1+3*L),]
  ci    <- output[rowname %in% (1+3*L+1):(1+4*L),]
  
  index <- rep(1:12, rep)
  temp.ptest <- temp.se <- temp.rej <- temp.ci  <- NULL
  for (i in 1:12) {
    temp.ptest <- rbind(temp.ptest, ptest[index == i,])
    temp.se    <- rbind(temp.se,    se[index == i,])
    temp.rej   <- rbind(temp.rej,   rej[index == i,])
    temp.ci    <- rbind(temp.ci,    ci[index == i,])
  }
  
  Bsq <- VAR <- MSE <- BtV <- CR <- IL <- matrix(NA, nrow=neval*4, ncol=ncol)
  true <- rep(tval, 4)
  for (i in 1:12) {
    Bsq[i,] <- (colMeans(temp.ptest[(1+(i-1)*rep): (i*rep),]) - true[i])^2
    VAR[i,] <- colVars(temp.ptest[(1+(i-1)*rep): (i*rep),])
    MSE[i,] <- (colMeans(temp.ptest[(1+(i-1)*rep): (i*rep),] - true[i])^2)
    BtV[i,] <- Bsq[i,] / VAR[i,]
    CR[i,]  <- (1 - colSums(temp.rej[(1+(i-1)*rep)  : (i*rep),]) / rep) * 100
    IL[i,]  <- colMeans(temp.ci[(1+(i-1)*rep)  : (i*rep),])
  }
  
  # columns: bias^2, variance, MSE, bias^2/variance, coverage rate, interval length
  # rows: classical estimator and three bias-corrected estimators
  table <- matrix(NA, nrow=ncol*4, ncol=6*neval)
  for (l in 1:4) {
    for (i in 1:neval) {
      table[(((l-1)*ncol+1):(l*ncol)),(((i-1)*6+1):(i*6))] <- cbind(Bsq[neval*(l-1)+i,],
                                                                    VAR[neval*(l-1)+i,],
                                                                    MSE[neval*(l-1)+i,],
                                                                    BtV[neval*(l-1)+i,],
                                                                    CR[neval*(l-1)+i,],
                                                                    IL[neval*(l-1)+i,])
    }
  }
  
  optk  <- as.matrix(read.table(paste("post/raw/optkappa_", j, ".txt", sep = ""), header = F))
  if (method == "wav") {
      ksel  <- round(colMeans(floor(log2(optk))), 1)
  } else {   
      ksel  <- round(colMeans(optk), 1)
  }

  subtable            <- table[,c(3, 5, 6, 9, 11, 12, 15, 17, 18)]
  subtable[,c(1,4,7)] <- sqrt(subtable[,c(1,4,7)])
  subtable            <- cbind(rep(c(range, ksel), 4), round(signif(subtable, 3), 3))
  if (method == "wav") subtable <- subtable[1:((length(range)+2)*3),]
  write.table(subtable, paste("post/pointwise/table_pointwise_", j, ".txt", sep = ""), row.names = F, col.names = F)
}


#######################
# Summary of kappa ####
#######################

# bs, evenly-spaced
sum.kappa <- NULL
for (j in 1:7) {
  optk  <- as.matrix(read.table(paste("post/raw/optkappa_", j, ".txt", sep = ""), header = F))
  sum.rot <- summary(optk[,1])
  sum.dpi <- summary(optk[,2])
  sum     <- cbind(rbind(sum.rot, sum.dpi), rbind(sd(optk[,1]), sd(optk[,2])))
  sum.kappa <- rbind(sum.kappa, sum)
  sum.kappa[, c(4, 7)] <- round(sum.kappa[, c(4, 7)], 2)
}
write.table(sum.kappa, paste("post/optk/sum_kappa_bs_es", ".txt", sep = ""), row.names = F, col.names = F)

#####################
# bs, quantile-spaced
sum.kappa <- NULL
for (j in 8:14) {
  optk  <- as.matrix(read.table(paste("post/raw/optkappa_", j, ".txt", sep = ""), header = F))
  sum.rot <- summary(optk[,1])
  sum.dpi <- summary(optk[,2])
  sum     <- cbind(rbind(sum.rot, sum.dpi), rbind(sd(optk[,1]), sd(optk[,2])))
  sum.kappa <- rbind(sum.kappa, sum)
  sum.kappa[, c(4, 7)] <- round(sum.kappa[, c(4, 7)], 2)
}
write.table(sum.kappa, paste("post/optk/sum_kappa_bs_qs", ".txt", sep = ""), row.names = F, col.names = F)

#################
# wav
sum.kappa <- NULL
for (j in 15:19) {
  optk  <- as.matrix(read.table(paste("post/raw/optkappa_", j, ".txt", sep = ""), header = F))
  optk <- floor(log2(optk))
  sum.rot <- summary(optk[,1])
  sum.dpi <- summary(optk[,2])
  sum     <- cbind(rbind(sum.rot, sum.dpi), rbind(sd(optk[,1]), sd(optk[,2])))
  sum.kappa <- rbind(sum.kappa, sum)
  sum.kappa[, c(4, 7)] <- round(sum.kappa[, c(4, 7)], 2)
}
write.table(sum.kappa, paste("post/optk/sum_kappa_wav", ".txt", sep = ""), row.names = F, col.names = F)
