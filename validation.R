# Compare with simulated data

# FDir values
fdir.red <- 0.862
fdir.nir <- 0.911
fdir.swir <- 0.966

# Soil reflectance
R.s.red <- 0.0825
R.s.nir <- 0.1363
R.s.swir <- 0.2139

# Leaf optics
rho.Ld.red <- 10.14/100
tau.Ld.red <- 5.26/100
rho.Ld.nir <- 45.25/100
tau.Ld.nir <- 49.13/100
rho.Ld.swir <- 31.66/100
tau.Ld.swir <- 41.63/100

# Other
LAI <- 2.2
theta.o <- 74.17
phi.o <- 81.9
ng <- 12

# Run model
source("run-disord.R")
mod.red <- run.disord(ng=ng, LAI=LAI, theta.o=theta.o, phi.o=phi.o,
                      R.s=R.s.red, rho.Ld=rho.Ld.red, tau.Ld=tau.Ld.red, fdir=fdir.red)
mod.nir <- run.disord(ng=ng, LAI=LAI, theta.o=theta.o, phi.o=phi.o,
                      R.s=R.s.nir, rho.Ld=rho.Ld.nir, tau.Ld=tau.Ld.nir, fdir=fdir.nir)
mod.swir <- run.disord(ng=ng, LAI=LAI, theta.o=theta.o, phi.o=phi.o,
                      R.s=R.s.swir, rho.Ld=rho.Ld.swir, tau.Ld=tau.Ld.swir, fdir=fdir.swir)
save(mod.red, mod.nir, mod.swir, file="validation.RData")

