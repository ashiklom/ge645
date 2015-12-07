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
hr.sz <- cbind(hr.sz.red, hr.sz.nir)
matplot(sz, hr.sz, type='l', lty=1, col=c("red", "blue"))

par(mfrow=c(1,2))
plot(sz, hr.sz.red, type='l', main="Figure 4.4: Red")
plot(sz, hr.sz.nir, type='l', main="Figure 4.4: NIR")

## Figure 4.5 Plot HDRF (RF.mat) as a function of phi.view
load("hdrf.RData")
ind <- 8 + (1:8)
ind2 <- c(5,6,7,8,1,2,3,4)
phi.v <- hdrf.run.a$RF.mat[ind, 2][ind2]
phi.v[1:4] <- phi.v[1:4] - 360
rf.a <- hdrf.run.a$RF.mat[ind, 3][ind2]
rf.b <- hdrf.run.b$RF.mat[ind, 3][ind2]
par(mfrow=c(1,2))
plot(phi.v, rf.a, type='l')
plot(phi.v, rf.b, type='l')
