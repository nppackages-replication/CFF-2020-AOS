#######################################################################
#### This file contains supporting functions used in simulation #######
#######################################################################
# Author: Matias D. Cattaneo, Max Farrell, Yingjie Feng
# Date:   05/19/2019

#####################################
########## Models ###################
#####################################
# Model 1
mu1 <- function(x) {
  x <- 2*x-1
  gx <- sin(pi*x/2) / (1 + 2 * x^2 *(sign(x)+1))
}

# Model 2
mu2 <- function(x) {
  x <- 2*x-1
  gx <- sin(3*pi*x/2) / (1 + 18* x^2 *(sign(x)+1))
}

# Model 3
mu3 <- function(x) {
  gx <- 2 * x - 1 + 5 * dnorm(20 * x - 10)
  return(gx)
}

# Model 4
mu4 <- function(x, y) {
  gx <- sin(5 * x) * sin(10 * y)
  return(gx)
}

# Model 5
mu5 <- function(x, y) {
  gx <- (1 - (4 * x - 2)^2)^2 * sin(5 * y) / 5
  return(gx)
}

# Model 6
mu6 <- function(x, y, z) {
  gx <- (1-(4*x-2)^2)^2 * (2*y-1) * (z-0.5)
  return(gx)
}

# Model 7
tau <- function(x) {
  gx <- (x-0.5)+8*(x-0.5)^2+6*(x-0.5)^3-30*(x-0.5)^4-30*(x-0.5)^5
  return(gx)
}

mu7 <- function(x, y, z) {
  gx <- tau(x) * tau(y) * tau(z)
  return(gx)
}

# List all funs
funlist <- list(mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, mu5=mu5,
                mu6=mu6, mu7=mu7)

################################################
############### DGP ############################
################################################
dgp <- function(n, model, d, sigma) {
  if (d == 2) {
    x <- matrix(runif(2*n), ncol = 2)
    y <- model(x[,1], x[,2]) + matrix(rnorm(n, 0, sigma), ncol = 1)
  } else if (d == 1) {
    x <- matrix(runif(n), ncol = 1)
    y <- model(x[,1]) + matrix(rnorm(n, 0, sigma), ncol = 1)
  } else if (d == 3) {
    x <- matrix(runif(3*n), ncol = 3)
    y <- model(x[,1], x[,2], x[,3]) + matrix(rnorm(n, 0, sigma), ncol = 1)
  }
  return(list(x = x, y = y))
}

################################################
############# Pointwise Inference ##############
################################################
compute <- function(data, eval, method, m, q, J, smooth, ktype, nknot, proj, d, model) {
  if (d == 3) {
    vce <- "hc1"
  } else {
    vce <- "hc2"
  }
  if (method != "wav") {
    Estimate <- lsprobust(data$y, data$x, eval, method = method, m = m, m.bc = q, proj = proj, vce=vce,
                          smooth = smooth, ktype = ktype, nknot = rep(nknot, d), bc = "bc1")$Estimate
    tau.cl  <- Estimate[, "tau.cl"]
    se.cl   <- Estimate[, "se.cl"]
    tau.bc1 <- Estimate[, "tau.bc"]
    se.bc1  <- Estimate[, "se.rb"]
    
    Estimate <- lsprobust(data$y, data$x, eval, method = method, m = m, m.bc = q, proj = proj, vce=vce,
                          smooth = smooth, ktype = ktype, nknot = rep(nknot, d), bc = "bc2")$Estimate
    tau.bc2 <- Estimate[, "tau.bc"]
    se.bc2  <- Estimate[, "se.rb"]
    
    Estimate <- lsprobust(data$y, data$x, eval, method = method, m = m, m.bc = q, proj = proj, vce=vce,
                           smooth = smooth, ktype = ktype, nknot = rep(nknot, d), bc = "bc3")$Estimate
    tau.bc3 <- Estimate[, "tau.bc"]
    se.bc3  <- Estimate[, "se.rb"]
    
  } else {
    Estimate <- lsprobust(data$y, data$x, eval, method = method, m = m, m.bc = q, J = rep(J, d), vce=vce,
                          bc = "bc1")$Estimate
    tau.cl  <- Estimate[, "tau.cl"]
    se.cl   <- Estimate[, "se.cl"]
    tau.bc1 <- Estimate[, "tau.bc"]
    se.bc1  <- Estimate[, "se.rb"]
    
    Estimate <- lsprobust(data$y, data$x, eval, method = method, m = m, m.bc = q, J = rep(J, d), vce=vce,
                          bc = "bc2")$Estimate
    tau.bc2 <- Estimate[, "tau.bc"]
    se.bc2  <- Estimate[, "se.rb"]
    
    tau.bc3 <- se.bc3 <- rep(NA, nrow(eval))
  }
  
  if (d == 3) mu.0 <- model(eval[,1], eval[,2], eval[,3])
  if (d == 2) mu.0 <- model(eval[,1], eval[,2])
  if (d == 1) mu.0 <- model(eval[,1])

  rej.cl  <- (((tau.cl  - mu.0) / se.cl  > qnorm(0.975)) | ((tau.cl  - mu.0) / se.cl  < qnorm(0.025))) * 1
  rej.bc1 <- (((tau.bc1 - mu.0) / se.bc1 > qnorm(0.975)) | ((tau.bc1 - mu.0) / se.bc1 < qnorm(0.025))) * 1
  rej.bc2 <- rej.bc3 <- rep(NA, nrow(eval))
  if (method != "pp") {
    rej.bc2 <- (((tau.bc2 - mu.0) / se.bc2 > qnorm(0.975)) | ((tau.bc2 - mu.0) / se.bc2 < qnorm(0.025))) * 1
    if (method == "bs") rej.bc3 <- (((tau.bc3 - mu.0) / se.bc3 > qnorm(0.975)) | ((tau.bc3 - mu.0) / se.bc3 < qnorm(0.025))) * 1
  }

  ci.cl   <- se.cl  * qnorm(0.975) * 2
  ci.bc1  <- se.bc1 * qnorm(0.975) * 2
  ci.bc2  <- ci.bc3 <- rep(NA, nrow(eval))

  if (method != "pp") {
    ci.bc2  <- se.bc2 * qnorm(0.975) * 2
    if (method == "bs") ci.bc3  <- se.bc3 * qnorm(0.975) * 2
  }

  result <- c(tau.cl, tau.bc1, tau.bc2, tau.bc3,
              se.cl,  se.bc1,  se.bc2,  se.bc3,
              rej.cl, rej.bc1, rej.bc2, rej.bc3,
              ci.cl,  ci.bc1,  ci.bc2,  ci.bc3)
  return(result)
}

sim <- function(i, n, d, eval, method, m, q, smooth, ktype, range, proj, model, sigma) {
  data <- dgp(n, model, d, sigma)
  if (d == 3) {
    vce <- "hc1"
  } else {
    vce <- "hc2"
  }
  optk <- lspkselect(data$y, data$x, m=m, method=method, vce=vce,
                     ktype=ktype, kselect="all", rotnorm=F)$ks[,c("k.rot", "k.dpi")]
  # avoid potential failure due to overfitting, not needed for current simuls
  if (d == 3 & method =="bs") {
    optk[optk>5] <- 5         
  }

  if (method == "wav") {
    wavk <- optk
    optk   <- pmax(floor(log2(optk)), ceiling(log2(2*(m+1))))
  }

  range <- c(range, optk)

  if (method == "bs" | method == "pp") {
    result <- sapply(range, function(nknot) compute(data, eval, method, m, q, J, smooth, ktype, nknot, proj, d, model))
  } else {
    result <- sapply(range, function(J) compute(data, eval, method, m, q, J, smooth, ktype, nknot, proj, d, model))
  }

  if (method == "wav") range[(length(range)-1):length(range)] <- wavk
  result <- rbind(range, result)
  return(result)
}



#####################################
##### Uniform Inference #############
#####################################
#####################################
unicov2 <- function(data, d, grid, method, m, q, J, smooth, ktype, nknot, proj, tval) {
  if (d == 3) {
    vce <- "hc1"
  } else {
    vce <- "hc2"
  }
  if (method != "wav") {
    # classical
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "none", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.cl"]
    se      <- fit$Estimate[, "se.cl"]
    cval    <- fit$sup.cval
    rej.pl.cl  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.cl <- mean(2*se*cval)
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "none", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.cl  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.cl <- mean(2*se*cval)
    
    # bc1
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc1", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    rej.pl.bc1  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.bc1 <- mean(2*se*cval)
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc1", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.bc1  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.bc1 <- mean(2*se*cval)
    
    # bc2
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc2", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    rej.pl.bc2  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.bc2 <- mean(2*se*cval)
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc2", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.bc2  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.bc2 <- mean(2*se*cval)
  
    # bc3
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc3", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    rej.pl.bc3  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.bc3 <- mean(2*se*cval)
  
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc3", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.bc3  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.bc3 <- mean(2*se*cval)
  
  } else {
    # classical
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "none", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.cl"]
    se      <- fit$Estimate[, "se.cl"]
    cval    <- fit$sup.cval
    rej.pl.cl  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.cl <- mean(2*se*cval)
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "none", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.cl  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.cl <- mean(2*se*cval)
    
    # bc1
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc1", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    rej.pl.bc1  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.bc1 <- mean(2*se*cval)
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc1", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.bc1  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.bc1 <- mean(2*se*cval)
    
    # bc2
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc2", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau     <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    rej.pl.bc2  <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.pl.bc2 <- mean(2*se*cval)
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc2", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    rej.wb.bc2   <- (((tau - tval) / se > cval)  | ((tau - tval) / se < -cval)) * 1
    width.wb.bc2 <- mean(2*se*cval)
    
    # bc3
    rej.pl.bc3 <- rep(NA, 100)
    rej.wb.bc3 <- rep(NA, 100)
    width.pl.bc3 <- NA
    width.wb.bc3 <- NA
  }

  result <- c(rej.pl.cl, rej.wb.cl,
              rej.pl.bc1, rej.wb.bc1,
              rej.pl.bc2, rej.wb.bc2,
              rej.pl.bc3, rej.wb.bc3,
              width.pl.cl, width.wb.cl,
              width.pl.bc1, width.wb.bc1,
              width.pl.bc2, width.wb.bc2,
              width.pl.bc3, width.wb.bc3)
}

sim.unicov2 <- function(i, n, d, model, grid, method, m, q, smooth, ktype, kimse, proj, tval, optk, sigma) {
  data  <- dgp(n, model, d, sigma)
  range <- c(kimse, optk[i,])
  if (method == "wav") range <- pmax(floor(log2(range)), ceiling(log2(2*(m+1))))
  kval  <- unique(range)

  if (method == "bs" | method == "pp") {
    result <- sapply(kval, function(nknot) unicov2(data, d, grid, method, m, q, J, smooth, ktype, nknot, proj, tval))
  } else {
    result <- sapply(kval, function(J) unicov2(data, d, grid, method, m, q, J, smooth, ktype, nknot, proj, tval))
  }

  output     <- matrix(NA, nrow=8*nrow(grid)+8, ncol=3)

  output[,1] <- result[,which(kval == range[1])]
  output[,2] <- result[,which(kval == range[2])]
  output[,3] <- result[,which(kval == range[3])]

  return(output)
}

banderror <- function(A, L, rep) {
  result <- sapply(1:rep, function(i) (colSums(A[((i-1)*L+1):(i*L),-1])>0)*1)
  result <- (1 - rowSums(result) / rep) * 100
  return(result)
}

pcoverror <- function(A, L, rep) {
  result <- sapply(1:L, function(i) colMeans(A[((i-1)*rep+1):(i*rep),]))
  return(t(result))
}

cpr <- function(A) {
  result <- colMeans(A <= 0.05)
  return(result)
}

###############################
##### Plot CB #################
###############################

plotband <- function(data, d, grid, method, m, q, J, smooth, ktype, nknot, proj) {
  if (d == 3) {
    vce <- "hc1"
  } else {
    vce <- "hc2"
  }
  if (method != "wav") {
    # classical
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "none", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau.cl  <- tau <- fit$Estimate[, "tau.cl"]
    se      <- fit$Estimate[, "se.cl"]
    cval    <- fit$sup.cval
    lb.pl.cl   <- tau - se*cval
    rb.pl.cl   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "none", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.cl   <- tau - se*cval
    rb.wb.cl   <- tau + se*cval
    
    # bc1
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc1", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    lb.pl.bc1   <- tau - se*cval
    rb.pl.bc1   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc1", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.bc1   <- tau - se*cval
    rb.wb.bc1   <- tau + se*cval
    
    # bc2
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc2", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    lb.pl.bc2   <- tau - se*cval
    rb.pl.bc2   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc2", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.bc2   <- tau - se*cval
    rb.wb.bc2   <- tau + se*cval
    
    # bc3
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc3", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    lb.pl.bc3   <- tau - se*cval
    rb.pl.bc3   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval=grid, method = method, m = m, m.bc = q,
                     proj = proj, smooth = smooth, ktype = ktype, nknot = rep(nknot, d), vce = vce,
                     bc = "bc3", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.bc3   <- tau - se*cval
    rb.wb.bc3   <- tau + se*cval
  } else {
    # classical
    fit <- lsprobust(data$y, data$x, eval = grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "none", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau.cl  <- tau <- fit$Estimate[, "tau.cl"]
    se      <- fit$Estimate[, "se.cl"]
    cval    <- fit$sup.cval
    lb.pl.cl   <- tau - se*cval
    rb.pl.cl   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval = grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "none", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.cl   <- tau - se*cval
    rb.wb.cl   <- tau + se*cval
    
    # bc1
    fit <- lsprobust(data$y, data$x, eval = grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc1", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    lb.pl.bc1   <- tau - se*cval
    rb.pl.bc1   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval = grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc1", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.bc1   <- tau - se*cval
    rb.wb.bc1   <- tau + se*cval
    
    # bc2
    fit <- lsprobust(data$y, data$x, eval = grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc2", uni.grid = grid, band = T, uni.method = "pl", B=1000, level=95)
    tau <- fit$Estimate[, "tau.bc"]
    se      <- fit$Estimate[, "se.rb"]
    cval    <- fit$sup.cval
    lb.pl.bc2   <- tau - se*cval
    rb.pl.bc2   <- tau + se*cval
    
    fit <- lsprobust(data$y, data$x, eval = grid, method = method, m = m, m.bc = q, vce = vce,
                     J = rep(J, d), bc = "bc2", uni.grid = grid, band = T, uni.method = "wb", B=1000, level=95)
    cval    <- fit$sup.cval
    lb.wb.bc2   <- tau - se*cval
    rb.wb.bc2   <- tau + se*cval
    
    lb.pl.bc3   <- rep(NA, nrows(grid))
    rb.pl.bc3   <- rep(NA, nrows(grid))
    lb.wb.bc3   <- rep(NA, nrows(grid))
    rb.wb.bc3   <- rep(NA, nrows(grid))
  }
  
  result <- c(tau.cl,
              lb.pl.cl, rb.pl.cl, lb.wb.cl, rb.wb.cl,
              lb.pl.bc1, rb.pl.bc1,lb.wb.bc1, rb.wb.bc1,
              lb.pl.bc2, rb.pl.bc2,lb.wb.bc2, rb.wb.bc2,
              lb.pl.bc3, rb.pl.bc3,lb.wb.bc3, rb.wb.bc3)
}
