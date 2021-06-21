##RSMs
yield <- function(xi1, xi2)
{
  xi1 <- 3*xi1 - 15
  xi2 <- xi2/50 - 13
  xi1 <- cos(0.5)*xi1 - sin(0.5)*xi2
  xi2 <- sin(0.5)*xi1 + cos(0.5)*xi2
  y <- exp(-xi1^2/80 - 0.5*(xi2 + 0.03*xi1^2 - 40*0.03)^2)
  return(100*y)
}

xi1 <- seq(1, 8, length=100)
xi2 <- seq(100, 1000, length=100)
g <- expand.grid(xi1, xi2)
y <- yield(g[,1], g[,2])
persp(xi1, xi2, matrix(y, ncol = length(xi2)),
      theta = 45,
      phi = 45,
      lwd=0.5,
      xlab = "xi1 : time",
      ylab = "xi2 : temperature",
      zlab = "yield",
      expand = 0.4)
par(mfrow=c(1,1))
cols <- heat.colors(128)
image(xi1, xi2, matrix(y, ncol = length(xi2)),
      col = cols,
      xlab =  "xi1 : time",
      ylab = "xi2 : temperature")
contour(xi1, xi2, matrix(y, ncol = length(xi2)),
        nlevels = 4,
        add = TRUE)

##Low-order Polynomials
first.order <- function(x1, x2)
{
  50 + 8*x1 + 3*x2
}

x1 <- x2 <- seq(-1, 1, length=100)
g <- expand.grid(x1, x2)
eta1 <- matrix(first.order(g[,1], g[,2]),
               ncol = length(x2))

##Aircraft Wing Weight

wingwt <- function(Sw=0.48,
                   Wfw=0.28,
                   A=0.38,
                   L=0.5,
                   q=0.62,
                   l=0.344,
                   Rtc=0.4,
                   Nz=0.37,
                   Wdg=0.38)
{
  ## put coded inputs back on natural scale
  Sw <- Sw*(200-150) + 150
  Wfw <- Wfw*(300 - 220) + 220
  A <- A*(10 - 6) + 6
  L <- (L*(10 - (-10)) - 10) * pi/180
  q <- q*(45 - 16) + 16
  l <- l*(1 - 0.5) + 0.5
  Rtc <- Rtc*(0.18 - 0.08) + 0.08
  Nz <- Nz*(6 - 2.5) + 2.5
  Wdg <- Wdg*(2500 - 1700) +1700
  
  ## calculation on natural scale
  W <- 0.036*Sw^0.758 * Wfw^0.0035 * (A/cos(L)^2)^0.6 * q^0.006
  W <- W * l^0.4 * (100*Rtc/cos(L))^(-0.3) * (Nz*Wdg)^(0.49)
  return(W)
}

## creating a grid
x <- seq(0, 1, length=100)
g <- expand.grid(x, x)

## using the grid to vary Nz and A
W.A.Nz <- wingwt(A=g[,1], Nz=g[,2])


## sets up a color palette that can be re-used from one experiment to the next
cs <- heat.colors(128)
bs <- seq(min(W.A.Nz), max(W.A.Nz), length=129)

image(x, x, matrix(W.A.Nz, ncol = length(x)),
      col = cs,
      breaks = bs,
      xlab = "A",
      ylab = "Nz")
contour(x, x, matrix(W.A.Nz,
                     ncol = length(x)),
        add = TRUE)

W.l.Wfw <- wingwt(l=g[,1],
                  Wfw=g[,2])


image(x, x, matrix(W.l.Wfw, ncol = length(x)),
      col = cs,
      breaks = bs,
      xlab = "l",
      ylab = "Wfw")

contour(x, x, matrix(W.l.Wfw, ncol = length(x)),
        add = TRUE)

#LHS
library(lhs)
n <- 1000
X <- data.frame(randomLHS(n, 9))
names(X) <- names(formals(wingwt))

plot(X[,1:2], pch=19, cex=0.5)
abline(h=c(0.6, 0.8), col=2, lwd=2)


Y <- wingwt(X[,1],
            X[,2],
            X[,3],
            X[,4],
            X[,5],
            X[,6],
            X[,7],
            X[,8],
            X[,9])

fit.lm <- lm(log(Y) ~ .^2, data = data.frame(Y,X))
fit.lmstep <- step(fit.lm, scope = formula(fit.lm),
                   direction = "backward",
                   k=log(length(Y)),
                   trace = 0)

coef(fit.lmstep)

#GP
library(laGP)

fit.gp <- newGPsep(X, Y, 2, 1e-6, dK=TRUE)
mle <- mleGPsep(fit.gp)

baseline <- matrix(rep(as.numeric(formals(wingwt)),
                       nrow(g)),
                   ncol = 9,
                   byrow = TRUE)
XX <- data.frame(baseline)
names(XX) <- names(X)

XX$A <- g[,1]
XX$Nz <- g[,2]

p <- predGPsep(fit.gp, XX, lite = TRUE)
image(x, x, matrix(p$mean, ncol = length(x)),
      col = cs,
      breaks = bs,
      xlab = "A",
      ylab = "Nz")
contour(x, x, matrix(p$mean, ncol = length(x)),
        add = TRUE)


## Rocket Boosters
lgbb1 <- read.table("lgbb/lgbb_original.txt", header = TRUE)
names(lgbb1)
nrow(lgbb1)

library(akima)
g <- interp(lgbb1$mach,
            lgbb1$alpha,
            lgbb1$lift,
            dupl="mean")

image(g, col = heat.colors(128),
      xlab = "mach",
      ylab = "alpha")
points(lgbb1$mach, lgbb1$alpha, cex=0.25, pch=18)

# Predicting Satellite drag
n <- 8
X <- data.frame(Umag=runif(n, 5500, 9500),
                Ts=runif(n, 100, 500),
                Ta=runif(n, 200, 2000),
                theta=runif(n, -pi, pi),
                phi=runif(n, -pi/2, pi/2),
                alphan=runif(n), sigmat=runif(n))
X <- rbind(X, X)
mesh <- "C:/Users/ebgor/tpm/tpm/Mesh_Files/GRACE_A0_B0_ascii_redone.stl"
moles <- c(0,0,0,0,1,0)

source("E:/Intergenerational Mobility/tpm/tpm/R/tpm.R")
system.time(y <- tpm(X, moles = moles, stl = mesh, verb = 0))

train <- read.csv("tpm/data/GRACE/CD_GRACE_1000_He.csv")
test <- read.csv("tpm/data/GRACE/CD_GRACE_100_He.csv")

r <- apply(rbind(train, test)[,1:7], 2, range)

### convert to coded inputs
X <- train[,1:7]
XX <- test[,1:7]
for(j in 1:ncol(X)) {
  X[,j] <- X[,j] - r[1,j]
  XX[,j] <- XX[,j] - r[1,j]
  X[,j] <- X[,j]/(r[2,j] - r[1,j])
  XX[,j] <- XX[,j]/(r[2,j] - r[1,j])
}

fit.gp <- newGPsep(X, train[,8], 2, 1e-6, dK=TRUE)
mle <- mleGPsep(fit.gp)
p <- predGPsep(fit.gp, XX, lite = TRUE)
rmspe <- sqrt(mean((100*(p$mean - test[,8])/test[,8])^2))
