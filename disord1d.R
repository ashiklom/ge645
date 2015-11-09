# Functions for discrete ordinates solution to RTE

# I.o.uncol.down {{{
# Downward uncollided direct solar radiation
I.o.uncol.down <- function(DeltaL, Nlayers, Gdir, I.o, mu.o){
    I.o.uc.d <- numeric(Nlayers)
    L1 <- 0
    for(k in 1:Nlayers){
        L2 <- k*DeltaL - (0.5*DeltaL)
        Prob <- exp(-1/abs(mu.o) * Gdir * (L2-L1))
        I.o.uc.d[k] <- I.o * Prob
    }
    return(I.o.uc.d)
}
# }}}

# soil.direct {{{
# Soil intensity
soil.direct <- function(DeltaL, Nlayers, Gdir, R.s){
    L1 <- 0
    L2 <- Nlayers * DeltaL
    F.o.uc.d.soil <- abs(mu.o) * I.o * exp(-1/abs(mu.o)) * Gdir * (L2-L1)
    I.o.uc.u.soil <- (R.s/pi) * F.o.uc.d.soil
    return(c("flux"=F.o.uc.d.soil, "intensity"=I.o.uc.u.soil))
}
# }}}

# I.o.uncol.up {{{
# Upward uncollided direct solar radiation
I.o.uncol.up <- function(DeltaL, Nlayers, Gdir, Gdif, I.o, mu.o,
                         ng, xg, R.s, I.o.uc.d.soil){
    I.o.uc.u <- array(0, c(Nlayers, ng, ng))
    L2 <- Nlayers * DeltaL
    for(i in (ng/2+1):ng){
        for(j in 1:ng){
            for(k in Nlayers:1){
                L1 <- k * DeltaL - (0.5 * DeltaL)
                Prob <- exp(-1/abs(xg[i]) * Gdif[j,i] * (L2-L1))
                I.o.uc.u[k,j,i] <- I.o.uc.u.soil * Prob
            }
        }
    }
    return(I.o.uc.u)
}
# }}}

# I.d.uncol.down {{{
# Downward uncollided diffuse sky radiation
I.d.uncol.down <- function(DeltaL, Nlayers, Gdif, I.d, ng, xg, I.d.uc.d){
    I.d.uc.d <- array(0, c(Nlayers, ng, ng))
    L1 <- 0
    for(i in 1:(ng/2)){
        for(j in 1:ng){
            for(k in 1:Nlayers){
                L2 <- k * DeltaL - (0.5*DeltaL)
                Prob <- exp(-1/abs(xg[i]) * Gdif[j,i] * (L2-L1))
                I.d.uc.d[k,j,i] <- I.d * Prob
            }
        }
    }
    return(I.d.uc.d)
}
# }}}

# soil.diffuse {{{
soil.diffuse <- function(DeltaL, Nlayers, ng, xg, wg, R.s){
    L1 <- 0
    L2 <- Nlayers * DeltaL
    sum1 <- 0
    for(i in 1:(ng/2)){
        upperlimit <- 2*pi
        lowerlimit <- 0
        conv <- (upperlimit-lowerlimit)/2
        sum2 <- 0
        for(j in 1:ng){
            Prob <- exp(-1/abs(xg[i]) * Gdif[j,i] * (L2-L1))
            I.d.uc.u.soil <- I.d * Prob
            sum2 <- sum2 + wg[j] * I.d.uc.u.soil
        }
        sum1 <- sum1 + wg[i] * abs(xg[j]) * sum2 * conv
    }
    F.d.uc.d.soil <- sum1
    I.d.uc.u.soil <- F.d.uc.d.soil * (R.s/pi)
    return(c("flux"=F.d.uc.d.soil, "intensity"=I.d.uc.u.soil))
}
# }}}

# I.d.uncol.up {{{
# Upward uncollided diffuse sky radiation
I.d.uncol.up <- function(DeltaL, Nlayers, Gdif, I.d, ng, xg, wg, R.s,
                         I.d.uc.u.soil){
    I.d.uc.u <- array(0, c(Nlayers, ng, ng))
    for(i in (ng/2+1):ng){
        for(j in 1:ng){
            L2 <- Nlayers * DeltaL
            for(k in Nlayers:1){
                L1 <- k * DeltaL - (0.5*DeltaL)
                Prob <- exp(-1/abs(xg[i]) * Gdif[j,i] * (L2-L1))
                I.d.uc.u[k,j,i] <- I.d.uc.u.soil * Prob
            }
        }
    }
    return(I.d.uc.u)
}
# }}}

# FCS {{{
# Evaluate first collision source
FCS <- function(Nlayers, ng, xg, wg, Gamma.d.dir, Gamma.d.dif,
         I.o.uc.d, I.o.uc.u, I.d.uc.d, I.d.uc.u){
    Q <- array(0, c(Nlayers, ng, ng))
    for(i in 1:ng){
        for(j in 1:ng){
            for(k in 1:Nlayers){
                Q[k,j,i] <- 1/pi * Gamma.d.dir[j,i] * I.o.uc.d[k]
            }
        }
    }
    for(i in 1:ng){
        for(j in 1:ng){
            for(k in 1:Nlayers){
                upperlimit1 <- 1
                lowerlimit1 <- -1
                conv11 <- (upperlimit1-lowerlimit1)/2
                sum1 <- 0
                for(n in 1:ng){
                    upperlimit2 <- 2 * pi
                    lowerlimit2 <- 0
                    conv21 <- (upperlimit2 - lowerlimit2)/2
                    sum2 <- 0
                    for(m in 1:ng){
                        sum2 <- sum2 + wg[m] * 1/pi * Gamma.d.dif[m,n,j,i] *
                            (I.o.uc.u[k,m,n] + I.d.uc.d[k,m,n] + I.d.uc.u[k,m,n])
                    }
                    sum1 <- sum1 * wg[n]*sum2*conv21
                }
                Q[k,j,i] <- Q[k,j,i] + sum1*conv11
            }
        }
    }
    return(Q)
}
# }}}

# SWEEP_DOWN {{{
SWEEP_DOWN <- function(Nlayers, ng, xg, wg, Gdif, DeltaL, JJ, R.s){
    Ic <- array(0, c(Nlayers+1, ng, ng))
    for(i in 1:(ng/2)){
        for(j in 1:ng){
            fij <- xg[i]/DeltaL - 0.5*Gdif[j,i]
            aij <- (0.5*Gdif[j,i] + xg[i]/DeltaL) / fij
            bij <- 1/fij
            for(k in 1:Nlayers){
                Ic[k+1,j,i] <- aij * Ic[k,j,i] - bij*JJ[k,j,i]
            }
        }
    }

# Evaluate flux density incident on ground
    sum1 <- 0
    for(i in 1:(ng/2)){
        upperlimit <- 2*pi
        lowerlimit <- 0
        conv <- (upperlimit - lowerlimit)/2
        sum2 <- 0
        for(j in 1:ng){
            sum2 <- sum2 + wg[j] * Ic[Nlayers+1,j,i]
        }
        sum1 <- sum1 + wg[i] * abs(xg[i])*sum2*conv
    }
    Fc.soil <- sum1

# Evaluate Ic upward at the ground
    for(i in (ng/2+1):ng){
        for(j in 1:ng){
            Ic[Nlayers+1,j,i] <- Fc.soil * R.s / pi
        }
    }
    return(Ic)
}
# }}}

# SWEEP_UP {{{
SWEEP_UP <- function(Nlayers, ng, xg, wg, Gdif, DeltaL, JJ, Ic){
    for(i in (ng/2+1):ng){
        for(j in 1:ng){
            fij <- xg[i]/DeltaL + 0.5*Gdif[j,i]
            cij <- (xg[i]/DeltaL - 0.5*Gdif[j,i]) / fij
            dij <- 1 / fij
            for(k in Nlayers:1){
                Ic[k,j,i] <- cij*Ic[k+1,j,i] + dij*JJ[k,j,i]
            }
        }
    }
    return(Ic)
}

# }}}

# check.convergence {{{
check.convergence <- function(Ic, Ic.old, epsilon){
    convergence <- TRUE
    for(i in (ng/2+1):ng){
        for(j in 1:ng){
            if(abs(Ic[1,j,i] - Ic.old[j,i]) > epsilon){
                convergence <- FALSE
                break
            }
        }
    }
    return(convergence)
}
# }}}

# MULTI_COLL_S {{{
MULTI_COLL_S <- function(Nlayers, ng, xg, wg, Gamma.d.dif, Ic){
    for(i in 1:ng){
        for(j in 1:ng){
            for(k in 1:Nlayers){
                upperlimit1 <- 1
                lowerlimit1 <- -1
                conv11 <- (upperlimit1 - lowerlimit1)/2
                sum1 <- 0
                for(n in 1:ng){
                    upperlimit2 <- 2*pi
                    lowerlimit2 <- 0
                    conv21 <- (upperlimit2 - lowerlimit2)/2
                    sum2 <- 0
                    for(m in 1:ng){
                        Ic.cell.center <- 0.5 * (Ic[k,m,n] + Ic[k+1,m,n])
                        sum2 <- sum2 + wg[m]*(1/pi) * Gamma.d.dif[m,n,j,i] *
                            Ic.cell.center
                    }
                    sum1 <- sum1 + wg[n] * sum2 * conv21
                }
                S[k,j,i] <- sum1 * conv11
            }
        }
    }
    return(S)
}
# }}}

# ENERGY_BAL {{{
# Compute the energy balance
ENERGY_BAL <- function(Nlayers, gq, mu.o, Q, S, DeltaL, R.s, rho.Ld, tau.Ld,
           Gdir, Gdif, F.o.uc.d.soil, F.d.uc.d.soil,
           I.o.uc.d, I.d.uc.d, I.o.uc.u, I.d.uc.u, Ic){
    upperlimit <- 2*pi
    lowerlimit <- 0
    conv <- (upperlimit-lowerlimit)/2
    HT.uc <- F.o.uc.d.soil + F.d.uc.d.soil
    sum1 <- 0
    for(i in 1:(ng/2)){
        sum2 <- 0
        for(j in 1:ng){
            sum2 <- sum2 + wg[j]*Ic[Nlayers+1,j,i]
        }
        sum1 <- sum1 + wg[i] * abs(xg[i]) * sum2 * conv
    }
    HT.c <- sum1

# Uncollided hemispherical reflectance
    L1 <- 0
    L2 <- Nlayers * DeltaL
    sum1 <- 0
    for(i in (ng/2+1):ng){
        sum2 <- 0
        for(j in 1:ng){
            Prob <- exp(-1/abs(xg[i]) * Gdif[j,i] * (L2-L1))
            I.uc.u <- (F.o.uc.d.soil + F.d.uc.d.soil) * R.s/pi * Prob
            sum2 <- sum2 + wg[j] * I.uc.u
        }
        sum1 <- sum1 + wg[i] * abs(xg[i]) * sum2 * conv
    }
    HR.uc <- sum1

    sum1 <- 0
    for(i in (ng/2+1):ng){
        sum2 <- 0
        for(j in 1:ng){
            sum2 <- sum2 + wg[j] * Ic[1,j,i]
        }
        sum1 <- sum1 + wg[i] * abs(xg[i]) * sum2 * conv
    }
    HR.c <- sum1

# Evaluate canopy absorption from I.o.uc.d
    sum1 <- 0
    for(k in 1:Nlayers){
        sum1 <- sum1 + I.o.uc.d[k]*Gdir
    }
    AB.o.uc.d <- sum1 * (1 - (rho.Ld + tau.Ld)) * DeltaL

# Evaluate canopy absorption from I.o.uc.u
    sum1 <- 0
    for(k in 1:Nlayers){
        for(i in (ng/2+1):ng){
            sum2 <- 0
            for(j in 1:ng){
                sum2 <- sum2 + wg[j]*I.o.uc.u[k,j,i] * Gdif[j,i]
            }
            sum1 <- sum1 + wg[i] * sum2 * conv
        }
    }
    AB.o.uc.u <- sum1 * (1 - (rho.Ld + tau.Ld)) * DeltaL

# Evaluate canopy absorption from I.d.uc.d
    sum <- 1
    for(k in 1:Nlayers){
        for(i in 1:(ng/2)){
            sum2 <- 0
            for(j in 1:ng){
                sum2 <- sum2 + wg[j] * I.d.uc.d[k,j,i] * Gdif[j,i]
            }
            sum1 <- sum1 + wg[i] * sum2 * conv
        }
    }
    AB.d.uc.d <- sum1 * (1 - (rho.Ld + tau.Ld)) * DeltaL

# Evaluate canopy absorption from I.d.uc.u
    sum1 <- 0
    for(k in 1:Nlayers){
        for(i in 1:(ng/2)){
            sum2 <- 0
            for(j in 1:ng){
                sum2 <- sum2 + wg[j] * I.d.uc.d[k,j,i] * Gdif[j,i]
            }
            sum1 <- sum1 + wg[i] * sum2 * conv
        }
    }
    AB.d.uc.u <- sum1 * (1 - (rho.Ld + tau.Ld)) * DeltaL

# Evaluate canopy absorption from I.c
    sum1 <- 0
    for(k in 1:Nlayers){
        for(i in 1:ng){
            sum2 <- 0
            for(j in 1:ng){
                sum2 <- sum2 + wg[j] * Ic[k,j,i] * Gdif[j,i]
            }
            sum1 <- sum1 + wg[i] * sum2 * conv
        }
    }
    AB.c <- sum1 * (1 - (rho.Ld + tau.Ld)) * DeltaL

    AB.uc <- AB.o.uc.d + AB.o.uc.u + AB.d.uc.d + AB.d.uc.u
    
    return(list(HR.uc=HR.uc, HR.c=HR.c, HT.uc=HT.uc, HT.c=HT.c,
                AB.uc=AB.uc, AB.c=AB.c))
}
# }}}
