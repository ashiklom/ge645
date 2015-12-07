source("run-disord.R")

# Figure 3
run1 <- run.disord(LAI=5,
                   theta.o=150,
                   phi.o=0,
                   fdir=0.7,
                   rho.Ld=0.475,
                   tau.Ld=0.45,
                   R.s=0.2)
save(run1, file="run1.RData")

# Figure 4.1 - 4.3
LAI <- seq(0.5, 5, 0.5) 
rho.red <- 0.075
tau.red <- 0.035
rs.red <- 0.125
rs.black <- 0
lai.list.red <- lapply(LAI, function(x) run.disord(LAI=x, rho.Ld=rho.red, tau.Ld=tau.red, R.s=rs.red,
                                                   theta.o=150, phi.o=0, fdir=0.7))
lai.list.red.bs <- lapply(LAI, function(x) run.disord(LAI=x, rho.Ld=rho.red, tau.Ld=tau.red, R.s=rs.black,
                                                   theta.o=150, phi.o=0, fdir=0.7))
rho.nir <- 0.475
tau.nir <- 0.45
rs.nir <- 0.2
lai.list.nir <- lapply(LAI, function(x) run.disord(LAI=x, rho.Ld=rho.nir, tau.Ld=tau.nir, R.s=rs.nir,
                                                   theta.o=150, phi.o=0, fdir=0.7))
lai.list.nir.bs <- lapply(LAI, function(x) run.disord(LAI=x, rho.Ld=rho.nir, tau.Ld=tau.nir, R.s=rs.black,
                                                   theta.o=150, phi.o=0, fdir=0.7))
save(LAI, lai.list.red, lai.list.red.bs, lai.list.nir, lai.list.nir.bs, file="lai.results.RData")

## Figure 4.4
sz <- c(seq(0, 80, length.out=8), 89)
sz.red <- lapply(sz, function(x) run.disord(theta.o=x, rho.Ld=rho.red, tau.Ld=tau.red, phi.o=0, fdir=0.7, LAI=3))
sz.nir <- lapply(sz, function(x) run.disord(theta.o=x, rho.Ld=rho.nir, tau.Ld=tau.nir, phi.o=0, fdir=0.7, LAI=3))
save(sz, sz.red, sz.nir, file="sz.results.RData")

## Figure 4.5 
rho.a <- 0.7
tau.a <- 0.225
hdrf.run.a <- run.disord(LAI=3, rho.Ld=rho.a, tau.Ld=tau.a, theta.o=15, phi.o=0, fdir=0.7, R.s=0.2)
rho.b <- 0.225
tau.b <- 0.7
hdrf.run.b <- run.disord(LAI=3, rho.Ld=rho.b, tau.Ld=tau.b, theta.o=15, phi.o=0, fdir=0.7, R.s=0.2)
save(hdrf.run.a, hdrf.run.b, file="hdrf.RData")
