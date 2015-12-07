#' ---
#' title: GE 645 RTM
#' author: Alexey Shiklomanov
#' output_format: pdf_document
#' ---

#' Check for equality
check.equal <- function(a, b, tol=1e-1)
    stopifnot(abs(a - b) < tol)

# gauss.quad {{{
#' Get Gauss quadrature of order `ng`. Return as list.
gauss.quad <- function(ng){
    xx <- c(-0.861136312,-0.339981044,-0.9324695,-0.6612094, 
            -0.2386192,-0.960289856,-0.796666477,-0.525532410,
            -0.183434642, -0.973906529,-0.865063367,-0.679409568,
            -0.433395394, -0.148874339, -0.981560634,-0.904117256,
            -0.769902674,-0.587317954, -0.367831499,-0.125233409)
    ww <- c(0.347854845, 0.652145155,0.1713245, 0.3607616,
            0.4679139, 0.101228536, 0.222381034, 0.313706646,
            0.362683783, 0.066671344,0.149451349,0.219086363,
            0.269266719, 0.295524225, 0.047175336,0.106939326,
            0.160078329,0.203167427, 0.233492537,0.249147046)
    ishift <- c(0,0,2,5,9,14)
    ng2 <- ng/2
    ng2s1 <- 1:ng2
    xg1 <- xx[ng2s1 + ishift[ng2]]
    wg1 <- ww[ng2s1 + ishift[ng2]]
    xg <- c(xg1, -rev(xg1))
    wg <- c(wg1, rev(wg1))
    out <- cbind(xg, wg)
    colnames(out) <- c("ordinates", "weights")
    return(out)
}
# }}}

# Check quad {{{
#' Check that sum(weights) = 2 and integral = 0.5.
check.quad <- function(gq){
    ng <- nrow(gq)
    s <- sum(gq[,"weights"])
    print(sprintf("Qweight is %f", s))
    check.equal(s, 2)
    p <- apply(gq, 1, prod)
    s2 <- sum(p[-(ng/2):0])
    print(sprintf("Qord = %f", s2))
    check.equal(s2, 1/2)
}
# }}}

# Leaf normal PDF {{{
#' Leaf normal PDF function (generic)
leaf.normal.pdf <- function(gq, fun){
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
    upperlimit <- pi/2
    lowerlimit <- 0
    conv1 <- (upperlimit - lowerlimit)/2
    conv2 <- (upperlimit + lowerlimit)/2
    neword <- conv1 * xg + conv2
    gL <- fun(neword)
    s <- sum(gL * wg) * conv1
    print(sprintf("Leaf normal PDF = %f", s))
    check.equal(s, 1)
    return(gL)
}
# }}}

# Some leaf PDFs {{{
planophile <- function(x) 2/pi * (1 + cos(2*x))
erectophile <- function(x) 2/pi * (1 - cos(2*x))
plagiophile <- function(x) 2/pi * (1 - cos(4*x))
extremophile <- function(x) 2/pi * (1 + cos(4*x))
uniform <- function(x) rep(2/pi, length(x))
spherical <- function(x) sin(x)
fun.list <- list(planophile, erectophile, plagiophile, extremophile, 
                 uniform, spherical)
# }}}

# G.dir.function {{{
#' Calculate G_FUNCTION given a direction of photon travel OMEGA^PRIME 
#' (muprime, phiprime) and the leaf normal PDF
G.dir.function <- function(gq, gL, hL, muprime, phiprime){
# Get quadrature
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[, "weights"]
# Declare parameters
    mu.tp <- muprime
    sin.tp <- sqrt(1 - muprime^2)
# Integrate over 0 to pi/2 for theta (vertical angles)
    upperlimit.tL <- pi/2
    lowerlimit.tL <- 0
    conv1.tL <- (upperlimit.tL - lowerlimit.tL) / 2
    conv2.tL <- (upperlimit.tL + lowerlimit.tL) / 2
# Integrate over 0 to 2pi for phi (azimuth)
    upperlimit.pL <- 2*pi
    lowerlimit.pL <- 0
    conv1.pL <- (upperlimit.pL - lowerlimit.pL) / 2
    conv2.pL <- (upperlimit.pL + lowerlimit.pL) / 2
# Integral over theta_L
    sum.tL <- 0
    for(i in 1:ng){
        neword.tL <- conv1.tL*xg[i] + conv2.tL
        mu.tL <- cos(neword.tL)
        sin.tL <- sin(neword.tL)
# Integral over phi_L
        sum.pL <- 0
        for(j in 1:ng){
            neword.pL <- conv1.pL*xg[j] + conv2.pL
            dotproduct <- abs(mu.tL*mu.tp + sin.tL*sin.tp * 
                              cos(neword.pL - phiprime))
            sum.pL <- sum.pL + wg[j]*hL[j] / (2*pi) * dotproduct
        }
        sum.pL <- sum.pL * conv1.pL
        sum.tL <- sum.tL + wg[i]*gL[i]*sum.pL
    }
    sum.tL <- sum.tL * conv1.tL
    Gdir <- sum.tL
    return(Gdir)
}
# }}}

# G.dif.function {{{
#' Evaluate G function in all quadrature directions and check for normalization
G.dif.function <- function(gq, gL, hL){
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
    upperlimit.pp <- 2*pi
    lowerlimit.pp <- 0
    conv1.pp <- (upperlimit.pp - lowerlimit.pp)/2
    conv2.pp <- (upperlimit.pp + lowerlimit.pp)/2
# Gdif matrix, direction by direction
    Gdif <- matrix(NA, ng, ng)
    for(i in 1:ng){
        muprime <- xg[i]
        for(j in 1:ng){
            phiprime <- conv1.pp*xg[j] + conv2.pp
            Gdif[j,i] <- G.dir.function(gq, gL, hL, muprime, phiprime)
        }
    }
# Check for normalization
    sum.tp <- 0
    for(i in (ng/2+1):ng){
        sum.pp <- 0
        for(j in 1:ng){
            sum.pp <- sum.pp + wg[j]*Gdif[j,i]
        }
        sum.pp <- sum.pp * conv1.pp
        sum.tp <- sum.tp + wg[i] * sum.pp
    }
    sum.tp <- sum.tp / (2*pi)
    print(sum.tp)
    check.equal(sum.tp, 0.5)
# Return
    return(Gdif)
}
# }}}

# Gamma.d.function {{{
Gamma.d.function <- function(gq, gL, hL, rho.Ld, tau.Ld,
                             muprime, phiprime, mu, phi){
# Extract quadrature
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
# Set parameters
    mu.t <- mu
    mu.tp <- muprime
    sin.t <- sqrt(1 - mu^2)
    sin.tp <- sqrt(1 - muprime^2)
# Limits for theta (zenith angle)
    upperlimit.tL <- pi/2
    lowerlimit.tL <- 0
    conv1.tL <- (upperlimit.tL - lowerlimit.tL)/2
    conv2.tL <- (upperlimit.tL + lowerlimit.tL)/2
# Limits for phi (azimuth angle)
    upperlimit.pL <- 2*pi
    lowerlimit.pL <- 0
    conv1.pL <- (upperlimit.pL - lowerlimit.pL)/2
    conv2.pL <- (upperlimit.pL + lowerlimit.pL)/2
# Integrate over theta_L
    sum.tL <- 0
    for(i in 1:ng){
        neword.tL <- conv1.tL*xg[i] + conv2.tL
        mu.tL <- cos(neword.tL)
        sin.tL <- sin(neword.tL)
# Integrate over phi_L
        sum.pL <- 0
        for(j in 1:ng){
            neword.pL <- conv1.pL*xg[j] + conv2.pL
            dotproduct1 <- mu.tL*mu.tp +
                sin.tL*sin.tp*cos(neword.pL - phiprime)
            dotproduct2 <- mu.tL*mu.t +
                sin.tL*sin.t*cos(neword.pL - phi)
            dotproduct <- dotproduct1*dotproduct2
            if (dotproduct < 0){
                sum.pL <- sum.pL + rho.Ld*wg[j]*hL[j]/(2*pi)*abs(dotproduct)
            } else {
                sum.pL <- sum.pL + tau.Ld*wg[j]*hL[j]/(2*pi)*abs(dotproduct)
            }
        }
        sum.pL <- sum.pL*conv1.pL
        sum.tL <- sum.tL + wg[i]*gL[i]*sum.pL
    }
    sum.tL <- sum.tL*conv1.tL
    Gamma.d <- sum.tL
    return(Gamma.d)
}
# }}}

# Gamma.d.dir.function {{{
#' This routine evaluates the GAMMA.d function (muprime, phiprime -> mu, phi) 
#' where (mu, phi) are quadrature directions and checks for normalization.
Gamma.d.dir.function <- function(gq, gL, hL, rho.Ld, tau.Ld, muprime, phiprime){
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
# Conversion factors
    upperlimit.p <- 2*pi
    lowerlimit.p <- 0
    conv1.p <- (upperlimit.p - lowerlimit.p)/2
    conv2.p <- (upperlimit.p + lowerlimit.p)/2
# Gamma.d.dir matrix direction by direction
    Gamma.d.dir <- matrix(0, ng, ng)
    for(i in 1:ng){
        mu <- xg[i]
        for(j in 1:ng){
            phi <- conv1.p*xg[j] + conv2.p
            Gamma.d.dir[j,i] <- Gamma.d.function(gq, gL, hL, rho.Ld, tau.Ld,
                                                 muprime, phiprime, mu, phi)
        }
    }
    return(Gamma.d.dir)
}
# }}}

#' GOOD THROUGH HERE

# check.Gamma.d.dir {{{
check.Gamma.d.dir <- function(gq, rho.Ld, tau.Ld, Gdir, Gamma.d.dir){
# Get quadrature
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
# Conversion factors for phiprime
    upperlimit.pp <- 2*pi
    lowerlimit.pp <- 0
    conv1.pp <- (upperlimit.pp - lowerlimit.pp)/2
    conv2.pp <- (upperlimit.pp + lowerlimit.pp)/2
# Check for normalization
    sum.tp <- 0
    for(i in 1:ng){
        sum.pp <- 0
        for(j in 1:ng){
            sum.pp <- sum.pp + wg[j]*Gamma.d.dir[j,i]
        }
        sum.pp <- sum.pp*conv1.pp
        sum.tp <- sum.tp + wg[i]*sum.pp
    }
    sum.tp <- sum.tp/pi
    sum.tp <- sum.tp/(Gdir*(rho.Ld+tau.Ld))
    print(sprintf("Gamma.d : 1.0 ?= %f", sum.tp))
    check.equal(sum.tp, 1)
}
# }}}

# Gamma.d.dif.function {{{
Gamma.d.dif.function <- function(gq, gL, hL, rho.Ld, tau.Ld, Gdif){
# Extract quadrature
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
# Conversion factors
    upperlimit.pp <- 2*pi
    lowerlimit.pp <- 0
    conv1.pp <- (upperlimit.pp - lowerlimit.pp)/2
    conv2.pp <- (upperlimit.pp + lowerlimit.pp)/2
# Gamma.d.dif matrix
    Gamma.d.dif <- array(0, rep(12,4))
    for(i in 1:ng){
        muprime <- xg[i]
        for(j in 1:ng){
            phiprime <- conv1.pp*xg[j] + conv2.pp
            dummy <- Gamma.d.dir.function(gq, gL, hL, rho.Ld, tau.Ld,
                                                muprime, phiprime)
            check.Gamma.d.dir(gq, rho.Ld, tau.Ld, Gdif[j,i], dummy)
            #Gamma.d.dif[j,i,,] <- dummy
            for(m in 1:ng){
                for(n in 1:ng){
                    Gamma.d.dif[j,i,n,m] <- dummy[n,m]
                }
            }
        }
    }
    return(Gamma.d.dif)
}
# }}}
