source("xsections.R")

#' Check quadrature and integration
gq <- gauss.quad(12)
check.quad(gq)
example.integral <- gq.integrate(function(x) x, 0, 1)
print(sprintf("Test integral = %f", example.integral))
check.equal(example.integral, 1/2)

#' Check leaf normal PDFs
plot(1, 1, type='n', xlim=c(0, pi/2), ylim=c(0,1.4))
for(i in 1:length(fun.list)){
    f <- fun.list[[i]]
    lnp <- leaf.normal.pdf(gq, f)
    x <- seq(0, pi/2, length.out=12)
    lines(x, lnp, col=i)
}
legend("left", c("plan", "erec", "plag", "extr", "unif", "spher"),
       col=1:6, lty=1)

#' Plot G function (Figure 8)
gL <- leaf.normal.pdf(gq, spherical)
hL <- rep(1, 12)
#hL <- leaf.normal.pdf(gq, uniform)
a <- 1
b <- 0
muprime <- seq(a, b, length.out=20)
phiprime <- pi/2
plot(1, 1, type='n', xlim=c(a,b), ylim=c(0.2,1), 
     main="Figure 8", xlab="mu", ylab="Geometry (G) function")
color <- 1
for(f in fun.list){
    gL <- leaf.normal.pdf(gq, f)
    Gdir <- sapply(muprime, function(x) G.dir.function(gq, gL, hL, x, phiprime))
    lines(muprime, Gdir, col=color)
    color <- color+1
}
legend("topleft", c("plan", "erec", "plag", "extr", "unif", "spher"),
       col=1:6, lty=1)
Gdif <- G.dif.function(gq, gL, hL)

#' Plot Gamma function (Figure 9)
muprime <- 1
phiprime <- pi/2
mu <- seq(-1, 1, length.out=12)
phi <- pi/2
plot(1,1,type='n', ylim=c(0,0.4), xlim=c(-1,1),
     xlab = "Cosine of scattering angle",
     ylab = "Gamma function",
     main = "Figure 9")
color <- 1
tseq <- seq(0, 0.5, length.out=5)
for(tw in tseq){
    w <- 1
    tau.Ld <- tw*w
    rho.Ld <- (1 - tw)*w
    Gamma.d <- sapply(mu, function(x)
                      Gamma.d.function(gq, gL, hL, rho.Ld, tau.Ld,
                                       muprime, phiprime, x, phi))
    lines(mu, Gamma.d, col=color, ylim=c(1,3))
    color <- color + 1
}
legend("topright", as.character(tseq), lty=1, col=1:5)

#' Run Gamma.d.dir.function
muprime <- cos(170 * pi/180)
phiprime <- 0
rho.Ld <- 0.5
tau.Ld <- 0.5
Gdir <- G.dir.function(gq, gL, hL, muprime, phiprime)
Gamma.g.dir <- Gamma.d.dir.function(gq, gL, hL, rho.Ld, tau.Ld, muprime, phiprime)
check.Gamma.d.dir(gq, rho.Ld, tau.Ld, Gdir, Gamma.g.dir)

#' Calculate full Gamma.d.dif matrix
Gdif <- G.dif.function(gq, gL, hL)
Gamma.d.dif <- Gamma.d.dif.function(gq, gL, hL, rho.Ld, tau.Ld, Gdif)
matplot(Gamma.d.dif[1,,,1], type='l')

#' Plot of normalized, azimuthally dependent phase function (Figure 11)
plot(1,1,type='n', xlim=c(-1,1), ylim=c(0,2),
     xlab = "Cosine of scattering angle",
     ylab = "P function",
     main = "Figure 11")
