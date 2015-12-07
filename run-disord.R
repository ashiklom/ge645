# Script to run discrete ordinates code
source("xsections.R")
source("disord1d.R")

# Get cross sections
run.disord <- function(maxiter = 100,
                       ng = 8, 
                       Nlayers = 100,
                       Epsilon = 0.0001,
                       theta.o = 120,
                       phi.o = 0,
                       Ftot = 1,
                       fdir = 0.7,
                       LAI = 3,
                       rho.Ld = 0.45,
                       tau.Ld = 0.45,
                       R.s = 0.3){
theta.o <- theta.o * pi/180
phi.o <- phi.o * pi/180
mu.o <- cos(theta.o)
I.o <- Ftot * fdir/abs(mu.o)
I.d <- Ftot * (1-fdir)/pi
DeltaL <- LAI/Nlayers

gq <- gauss.quad(ng)
xg <- gq[,"ordinates"]
wg <- gq[, "weights"]
gL <- leaf.normal.pdf(gq, planophile)
hL <- rep(1, ng)
Gdir <- G.dir.function(gq, gL, hL, mu.o, phi.o)
Gdif <- G.dif.function(gq, gL, hL)
Gamma.d.dir <- Gamma.d.dir.function(gq, gL, hL, rho.Ld, tau.Ld, mu.o, phi.o)
Gamma.d.dif <- Gamma.d.dif.function(gq, gL, hL, rho.Ld, tau.Ld, Gdif)

# Evaluate uncollided direct solar radiation (I.o.uc.d, I.o.uc.u)
I.o.uc.d <- I.o.uncol.down(DeltaL, Nlayers, Gdir, I.o, mu.o)
soil.dir <- soil.direct(DeltaL, Nlayers, I.o, mu.o, Gdir, R.s)
F.o.uc.d.soil <- soil.dir["flux"]
I.o.uc.u.soil <- soil.dir["intensity"]
I.o.uc.u <- I.o.uncol.up(DeltaL, Nlayers, Gdir, Gdif, I.o, mu.o,
                         ng, xg, R.s, I.o.uc.u.soil)

# Evaluate uncollided diffuse sky radiation (I.d.uc.d, I.d.uc.u)
I.d.uc.d <- I.d.uncol.down(DeltaL, Nlayers, Gdif, I.d, ng, xg)
soil.dif <- soil.diffuse(DeltaL, Nlayers, ng, xg, wg, Gdif, I.d, R.s)
F.d.uc.d.soil <- soil.dif["flux"]
I.d.uc.u.soil <- soil.dif["intensity"]
I.d.uc.u <- I.d.uncol.up(DeltaL, Nlayers, Gdif, I.d, ng, xg, wg, R.s, I.d.uc.u.soil)

# Evaluate first collision source
Q <- FCS(Nlayers, ng, xg, wg, Gamma.d.dir, Gamma.d.dif,
         I.o.uc.d, I.o.uc.u, I.d.uc.d, I.d.uc.u)

# Iterate on Multiple-Collision Source S
S <- array(0, c(Nlayers, ng, ng))
S.store <- list()
Ic <- array(0, c(Nlayers+1, ng, ng))
for(ims in 1:maxiter){
    print(ims)
    S <- Q + S

# Sweep downwards plus handle the bottom boundary condition
    Ic <- SWEEP_DOWN(Nlayers, ng, xg, wg, Gdif, DeltaL, S, R.s, Ic)

# Sweep upwards and check for convergence
    Ic.old <- Ic[1,,]
    Ic <- SWEEP_UP(Nlayers, ng, xg, wg, Gdif, DeltaL, S, Ic) 
    convergence <- check.convergence(ng, Ic, Ic.old, Epsilon)

# Evaluate Multiple-Collision source
    if(convergence) break 
    S <- MULTI_COLL_S(Nlayers, ng, xg, wg, Gamma.d.dif, Ic)

# Store collision source for each iteration
    S.store[[ims]] <- S
}

# Do energy balance
eb.list <- ENERGY_BAL(Nlayers, ng, xg, wg, mu.o, Q, S, DeltaL, R.s, rho.Ld, tau.Ld,
           Gdir, Gdif, F.o.uc.d.soil, F.d.uc.d.soil,
           I.o.uc.d, I.d.uc.d, I.o.uc.u, I.d.uc.u, Ic)
## eb.list consists of: HR.uc, HR.c, HT.uc, HT.c, AB.uc, AB.c
eb.list[["HR"]] <- with(eb.list, HR.uc + HR.c)
eb.list[["HT"]] <- with(eb.list, HT.uc + HT.c)
eb.list[["AB"]] <- with(eb.list, AB.uc + AB.c)

eb.list[["USA"]] <- with(eb.list, (1-R.s) * HT.c)
eb.list[["CSA"]] <- with(eb.list, (1-R.s) * HT.c)
eb.list[["SA"]] <- with(eb.list, (1-R.s) * HT)
eb.list[["EB"]] <- with(eb.list, HR + AB + (1-R.s) * HT)

RF.mat <- matrix(NA, nrow = ng/2 * ng, ncol=3)
colnames(RF.mat) <- c("theta.v", "phi.v", "RF")
x <- 1
for(i in (ng/2+1):ng){
    theta.v <- xg[i] * 180/pi
    for(j in 1:ng){
        phi.v <- xg[j]*pi + pi
        phi.v <- phi.v * 180/pi
        RF <- (Ic[1,j,i] + I.o.uc.u[1,j,i] + I.d.uc.u[1,j,i]) * pi
        RF.mat[x,] <- c(theta.v, phi.v, RF)
        x <- x + 1
        #print(paste("Theta.view, Phi.view, Refl. factor:", 
                    #theta.v, phi.v, RF))
    }
}

out <- list(S.store = S.store, eb.list = eb.list, RF.mat = RF.mat)
return(out)
}


#with(eb.list, {
         #print(paste0("Uncollided hemispherical reflectance: ", HR.uc))
         #print(paste0("Collided Hemispherical reflectance: ", HR.c))
         #print(paste0("Hemispherical Reflectance: ", HR))
         #print(paste0("Uncollided Hemispherical Transmittance: ", HT.uc))
         #print(paste0("Collided Hemispherical Transmittance: ", HT.c))
         #print(paste0("Hemispherical Transmittance: ", HT))
         #print(paste0("Uncollided Canopy Absorbance: ", AB.uc))
         #print(paste0("Collided Canopy Absorbance: ", AB.c))
         #print(paste0("Canopy Absorbance: ", AB))
         #print(paste0("Uncollided Soil Absorbance: ", (1-R.s)*HT.uc))
         #print(paste0("Collided Soil Absorbance: ", (1-R.s)*HT.c))
         #print(paste0("Soil Absorbance: ", (1-R.s)*HT))
         #print(paste0("Energy balance (=1?): ", HR + AB + (1-R.s)*HT))
#})


