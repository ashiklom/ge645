#' ---
#' title: GE 645 RTM
#' author: Alexey Shiklomanov
#' output_format: pdf_document
#' ---

#' Check for equality
check.equal <- function(a, b, tol=1e-2)
    stopifnot(abs(a - b) < tol)

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
gq <- gauss.quad(12)

#' Check that sum(weights) = 2 and integral = 0.5.

tol <- 1e-2
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
check.quad(gq)

#' Calculate values at ordinates
gq.calc <- function(fun, lowerlimit, upperlimit){
    conv1 <- (upperlimit - lowerlimit)/2
    conv2 <- (upperlimit + lowerlimit)/2
    neword <- conv1*gq[,"ordinates"] + conv2
    y <- fun(neword)
    return(y)
}

#' Generic integration
gq.integrate <- function(fun, lowerlimit, upperlimit){
    conv1 <- (upperlimit - lowerlimit)/2
    conv2 <- (upperlimit + lowerlimit)/2
    neword <- conv1*gq[,"ordinates"] + conv2
    y <- fun(neword)
    s <- sum(y * gq[, "weights"]) * conv1
    return(s)
}

example.integral <- gq.integrate(function(x) x, 0, 1)
print(sprintf("Test integral = %f", example.integral))
check.equal(example.integral, 1/2)

#' Leaf normal PDF
leaf.normal.pdf <- function(gq, fun){
    ng <- nrow(gq)
    upperlimit <- pi/2
    lowerlimit <- 0
    gL <- gq.calc(fun, lowerlimit, upperlimit)
    s <- gq.integrate(fun, lowerlimit, upperlimit)
    print(sprintf("Leaf normal PDF = %f", s))
    check.equal(s, 1)
    return(gL)
}

#' Some functions:
planophile <- function(x) 2/pi * (1 + cos(2*x))
erectophile <- function(x) 2/pi * (1 - cos(2*x))
plagiophile <- function(x) 2/pi * (1 - cos(4*x))
extremophile <- function(x) 2/pi * (1 + cos(4*x))
uniform <- function(x) rep(2/pi, length(x))
spherical <- function(x) sin(x)

fun.list <- list(planophile, erectophile, plagiophile, extremophile, 
                 uniform, spherical)

plot(1, 1, type='n', xlim=c(0, pi/2), ylim=c(0,1.4))
for(i in 1:length(fun.list)){
    f <- fun.list[[i]]
    lnp <- leaf.normal.pdf(gq, f)
    x <- seq(0, pi/2, length.out=12)
    lines(x, lnp, col=i)
}
legend("left", c("plan", "erec", "plag", "extr", "unif", "spher"),
       col=1:6, lty=1)

#' G_function for a direction of photon travel
G_dir_function <- function(gq, gL, hL, muprime, phiprime){
    ng <- nrow(gq)
    mu.tp <- muprime
    sin.tp <- sqrt(1 - muprime^2)
# Integrate over 0 to pi/2 for theta (vertical angles)
    upperlimit.tL <- pi/2
    lowerlimit.tL <- 0
    conv1.tL <- (upperlimit.tL - lowerlimit.tL) / 2
    conv2.tL <- (upperlimit.tL + lowerlimit.tL) / 2
# Integrate over 0 to 2pi for phi (azimuth)
    upperlimit.pL <- 2 * pi
    lowerlimit.pL <- 0
    conv1.pL <- (upperlimit.pL - lowerlimit.tL) / 2
    conv2.pL <- (upperlimit.pL + lowerlimit.pL) / 2
    neword.tL <- conv1.tL * gq[,"ordinates"] + conv2.tL
    neword.pL <- conv1.pL * gq[,"ordinates"] + conv2.pL
    mu.tL <- cos(neword.tL)
    sin.tL <- sin(neword.tL)
# Perform integration in matrix form
    dotproduct <- abs(matrix(sin.tL, ncol=1) %*% cos(neword.pL - phiprime) *
        sin.tp + mu.tL*mu.tp)
    mm <- matrix(gq[,"weights"] * hL / (2*pi) * conv1.pL) %*%
        (gq[,"weights"] * gL * conv1.tL) * dotproduct
    Gdir <- sum(mm)
    return(Gdir)
}
gL <- leaf.normal.pdf(gq, spherical)
hL <- rep(1, nrow(gq))
muprime <- 0
phiprime <- 0
Gdir <- G_dir_function(gq, gL, hL, muprime, phiprime)
print(gd)


#' Evaluate G function in all quadrature directions and check for normalization

G_dif_function <- function(gq, gL, hL){
    upperlimit.pp <- 2*pi
    lowerlimit.pp <- 0
    conv1.pp <- (upperlimit.pp - lowerlimit.pp)/2
    conv2.pp <- (upperlimit.pp + lowerlimit.pp)/2
# Gdif matrix
    muprime <- gq[,"ordinates"]
    phiprime <- conv1.pp * gq[,"ordinates"] + conv2.pp
    prime <- expand.grid(muprime, phiprime)
    Gdif.raw <- apply(prime, 1, function(x) G_dir_function(gq, gL, hL, x[1], x[2]))
    Gdif <- matrix(Gdif.raw, 12, byrow=TRUE)
# Validate Gdif matrix -- integral over upper hemisphere should sum to 0.5
    ng2 <- (ng/2+1):ng
    check.pp <- Gdif[,ng2] * gq[,"weights"] * conv1.pp
    check.tp <- t(t(check.pp) * gq[ng2, "weights"])
    sum.tp <- sum(check.tp)/(2 * pi)
    check.equal(sum.tp, 0.5)
    return(Gdif)
}
Gdif <- G_dif_function(gq, gL, hL)

#' Calculate G_FUNCTION given a direction of photon travel OMEGA^PRIME and the 
#' leaf normal PDF
Gamma.d.function <- function(gq, gL, hL, rho.Ld, tau.Ld,
                             muprime, phiprime, mu, phi){
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
    mu.t <- mu
    mu.tp <- muprime
    sin.t <- sqrt(1 - mu^2)
    sin.tp <- sqrt(1 - muprime^2)
# Limits for theta (zenith angle)
    upperlimit.tL <- pi/2
    lowerlimit.tL <- 0
    conv1.tL <- (upperlimit.tL - lowerlimit.tL)/2
    conv2.tl <- (upperlimit.tL + lowerlimit.tL)/2
# Limits for phi (azimuth angle)
    upperlimit.tP <- 2*pi
    lowerlimit.tP <- 0
    conv1.tP <- (upperlimit.tP - lowerlimit.tP)/2
    conv2.tl <- (upperlimit.tP + lowerlimit.tP)/2
# Double integral in matrix form
    neword.tL <- conv1.tL * xg + conv2.tL
    neword.pL <- conv1.pL * xg + conv2.pL
    mu.tL <- cos(neword.tL)
    sin.tL <- sin(neword.tL)
    dotproduct1 <- matrix(sin.tL) %*% cos(neword.pL - phiprime) *
        sin.tp + mu.tL*mu.tp
    dotproduct2 <- matrix(sin.tL) %*% cos(neword.pL - phi) *
        sin.t + mu.tL*mu.t
# Code for checking dotproducts are calculated correctly
    #i <- 7
    #j <- 10
    #dp1 <- mu.tL[i] * mu.tp + sin.tL[i] * sin.tp * cos(neword.pL[j] - phiprime)
    #print(c(dp1, dotproduct1[i,j]))
    #dp2 <- mu.tL[i] * mu.t + sin.tL[i] * sin.t * cos(neword.pL[j] - phi)
    #print(c(dp2, dotproduct2[i,j]))
    dotproduct <- dotproduct1 * dotproduct2
    dp.neg <- dotproduct < 0
    dp.pos <- !dp.neg
    RT <- dp.neg * rho.Ld + dp.pos * tau.Ld
    pL.mat <- t(t(abs(dotproduct) * RT) * wg * hL) / (2*pi) * conv1.pL
    tL.mat <- pL.mat * wg * gL
    #pL.mat <- RT * wg * hL / (2*pi) * abs(dotproduct) * conv1.pL
    #tL.mat <- t(t(pL.mat) * wg * gL)
    Gamma.d <- sum(tL.mat) * conv1.tL
    return(Gamma.d)
}
rho.Ld <- 0.4
tau.Ld <- 0.0
mu <- seq(-1, 1, length.out=12)
phi <- 0
Gamma.d <- sapply(mu, function(x)
                Gamma.d.function(gq, gL, hL, rho.Ld, tau.Ld,
                                 muprime, phiprime, x, phi))
plot(mu, Gamma.d, type='l')

#' This routine evaluates the GAMMA.d function (muprime, phiprime -> mu, phi) 
#' where (mu, phi) are quadrature directions and checks for normalization.

Gamma.d.dir.function <- function(gq, gL, hL, rho.Ld, tau.Ld, muprime, phiprime){
    ng <- nrow(gq)
    xg <- gq[,"ordinates"]
    wg <- gq[,"weights"]
    upperlimit.p <- 2*pi
    lowerlimit.p <- 0
    conv1.p <- (upperlimit.p - lowerlimit.p)/2
    conv2.p <- (upperlimit.p + lowerlimit.p)/2
# Gamma.d.dir matrix direction by direction
    mu <- xg
    phi <- conv1.p * xg + conv2.p
    mu.phi <- expand.grid(mu, phi)
    Gamma.d.raw <- apply(prime, 1, 
                         function(x) Gamma.d.function(gq, gL, hL,
                                                      rho.Ld, tau.Ld,
                                                      muprime, phiprime,
                                                      x[1], x[2]))
    Gamma.d.dir <- matrix(Gamma.d.raw, ng, byrow=TRUE)
# Check -- should normalize to 1
    Gdir <- G_dir_function(gq, gL, hL, muprime, phiprime)
    check.pp <- Gamma.d.dir * wg * conv1.p
    check.tp <- t(t(check.pp) * wg)
    sum.tp <- sum(check.tp)/(pi * Gdir * (rho.Ld + tau.Ld))
    check.equal(sum.pp, 1)
    return(Gamma.d.dir)
}

Gamma.d.dif.function <- function(gq, gL, hL, rho.Ld, tau.Ldh
    
    
                                 

###


