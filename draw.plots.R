#' Run cross sections plots
source("test-xsections.R")

load("run1.RData")
S.store <- run1$S.store

# Figure 3 -- S at top and bottom as a function of iteration
S.plot <- sapply(S.store, "[", c(1,100), 1, 1)
matplot(t(S.plot), type='l', xlab="Iteration", ylab="Multiple scattering source")

load("lai.results.RData")

hr.red <- sapply(lai.list.red, function(x) x$eb.list$HR)
hr.red.bs <- sapply(lai.list.red.bs, function(x) x$eb.list$HR)
hr.nir <- sapply(lai.list.nir, function(x) x$eb.list$HR)
hr.nir.bs <- sapply(lai.list.nir.bs, function(x) x$eb.list$HR)
ht.red <- sapply(lai.list.red, function(x) x$eb.list$HT)
ht.nir <- sapply(lai.list.nir, function(x) x$eb.list$HT)

plot(LAI, hr.red, type='l', ylim=c(0.02,0.09)) ## 4.1
lines(LAI, hr.red.bs, lty=2)
legend("topright", c("Refl.", "Black"), lty=1:2, title="Soil")
title(main="Figure 4.1: Red", xlab="LAI", ylab="Bi-hemispherical reflectance")

plot(LAI, hr.nir, type='l', ylim=c(0.15,0.55)) ## 4.2
lines(LAI, hr.nir.bs, lty=2)
legend("bottomright", c("Refl.", "Black"), lty=1:2, title="Soil")
title(main="Figure 4.2: NIR", xlab="LAI", ylab="Bi-hemispherical reflectance")

## 4.3
plot(LAI, ht.red, type='l', col="red")
lines(LAI, ht.nir, col="blue")
legend("topright", c("Red", "NIR"), lty=1, col=c("red", "blue"))
title(main="Figure 4.3: Transmittance", xlab="LAI", ylab="Bi-hemispherical transmittance")


## Figure 4.4 Plot HR for red (a) and NIR (b) as a function of theta.o (sz)
load("sz.results.RData")
hr.sz.red <- sapply(sz.red, function(x) x$eb.list$HR)
hr.sz.nir <- sapply(sz.nir, function(x) x$eb.list$HR)
#hr.sz <- cbind(hr.sz.red, hr.sz.nir)
#matplot(sz, hr.sz, type='l', lty=1, col=c("red", "blue"))
par(mfrow=c(1,2))
plot(sz, hr.sz.red, type='l', main="Figure 4.4: Red")
plot(sz, hr.sz.nir, type='l', main="Figure 4.4: NIR")

## Figure 4.5 Plot HDRF (RF.mat) as a function of phi.view
load("hdrf.RData")

hdrf.plot <- function(rf.mat, ng, i1, i2){
# Backscatter when azimuth ~ 0
    nr <- nrow(rf.mat)
    theta.fw <- -rev(unique(rf.mat[,1]))
    ind.fw <- seq(i1, nr, ng)
    hdrf.fw <- rf.mat[ind.fw, 3]
# Forwardscatter when azimuth ~ 180
    theta.bw <- -theta.fw
    ind.bw <- seq(i2, nr, ng)
    hdrf.bw <- rf.mat[ind.bw,3]
# Combine into single matrix
    theta <- c(theta.bw, theta.fw)
    hdrf <- c(hdrf.bw, hdrf.fw)
    hdrf.mat <- cbind(theta, hdrf)
    hdrf.mat <- hdrf.mat[order(theta),]
    return(hdrf.mat)
}
hdrf.a <- hdrf.plot(hdrf.run.a$RF.mat, 8, 1, 5)
hdrf.b <- hdrf.plot(hdrf.run.b$RF.mat, 8, 1, 5)
par(mfrow=c(1,1))
plot(hdrf.a, type='l', col="blue", ylim=c(0.44, 0.6),
     main= "Figure 4.5")
lines(hdrf.b, type='l', col="red")
legend("topright", c("a", "b"), lty=1, col=c("red", "blue"))

# Validation plots
load("validation.RData")
dat <- read.csv("pairie_site_validation data.csv")
dat[,-1] <- dat[,-1] / 100      # Convert from percent reflectance to fraction

# Validation model inputs
hdrf.red <- hdrf.plot(mod.red$RF.mat, 12, 1, 7)
hdrf.red.2 <- hdrf.plot(mod.red$RF.mat, 12, 3, 9)
par(mfrow=c(1,2))
plot(red1 ~ view.zenith, data=dat, col="red", pch=2, ylim=c(0, 0.15))
points(hdrf.red, col="blue", pch=3)
plot(red1 ~ view.zenith, data=dat, col="red", pch=2, ylim=c(0, 0.15))
points(hdrf.red.2, col="blue", pch=3)

hdrf.nir <- hdrf.plot(mod.nir$RF.mat, 12, 1, 7)
hdrf.nir.2 <- hdrf.plot(mod.nir$RF.mat, 12, 3, 9)
par(mfrow=c(1,2))
plot(nir1 ~ view.zenith, data=dat, col="red", pch=2, ylim=c(0, 0.8))
points(hdrf.nir, col="blue", pch=3)
plot(nir1 ~ view.zenith, data=dat, col="red", pch=2, ylim=c(0, 0.8))
points(hdrf.nir.2, col="blue", pch=3)

hdrf.swir <- hdrf.plot(mod.swir$RF.mat, 12, 1, 7)
hdrf.swir.2 <- hdrf.plot(mod.swir$RF.mat, 12, 3, 9)
par(mfrow=c(1,2))
plot(swir1 ~ view.zenith, data=dat, col="red", pch=2, ylim=c(0, 0.7))
points(hdrf.swir, col="blue", pch=3)
plot(swir1 ~ view.zenith, data=dat, col="red", pch=2, ylim=c(0, 0.7))
points(hdrf.swir.2, col="blue", pch=3)
