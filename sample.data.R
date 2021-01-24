# simulate data like Saaerela
# msmsim('./results0.txt', from=1, to=2, nobs=2000, ncov=5, coeffa=0.3, coeffb=0.3, coeffc=-0.75, effect=-0.25, a0=-1.0, dobayes=TRUE, doboot=TRUE)
library(mvtnorm)
set.seed(1)

ncov = 5 # number of (time-varying) covariates
a0=-1
coeffa=0.3; coeffb=0.3; coeffc=-0.75; effect=-0.25
nobs=500
ntot = 1000 # total number of patients

x <- array(NA, c(2, 2, ntot, 3 * ncov)) # theoretical matrix/ 3 time periods
xobs <- matrix(NA, ntot, 3 * ncov) # observed matrix
y <- array(NA, c(2, 2, 2, ntot)) # counter-factual response
z <- matrix(NA, ntot, 3) # treatment
xdiffsd <- 0.25
a <- rep(coeffa, ncov)
b <- rep(coeffb, ncov)
s <- matrix(0.25, ncov, ncov)
diag(s) <- 1.0
su <- matrix(0.05, ncov, ncov)
diag(su) <- 0.1


# 1. intervention:

u <- rmvnorm(ntot, rep(0.0, ncov), su)
xobs[,1:ncov] <- x[1,1,,1:ncov] <- x[1,2,,1:ncov] <- x[2,1,,1:ncov] <- x[2,2,,1:ncov] <- rmvnorm(ntot, rep(0.0, ncov), s)
pz <- 1.0/(1.0+exp(-(a0 + xobs[,1:ncov] %*% a)))
z[,1] <- (runif(ntot) < pz)

# 2. intervention:

x[1,1,,(ncov+1):(2*ncov)] <- x[1,2,,(ncov+1):(2*ncov)] <- rmvnorm(ntot, rep(0.0, ncov), (xdiffsd^2) * s)
x[1,1,,(ncov+1):(2*ncov)] <- x[1,1,,(ncov+1):(2*ncov)] + x[1,1,,1:ncov] + u
x[1,2,,(ncov+1):(2*ncov)] <- x[1,2,,(ncov+1):(2*ncov)] + x[1,2,,1:ncov] + u

x[2,1,,(ncov+1):(2*ncov)] <- x[2,2,,(ncov+1):(2*ncov)] <- rmvnorm(ntot, rep(coeffc, ncov), (xdiffsd^2) * s)
x[2,1,,(ncov+1):(2*ncov)] <- x[2,1,,(ncov+1):(2*ncov)] + x[2,1,,1:ncov] + u
x[2,2,,(ncov+1):(2*ncov)] <- x[2,2,,(ncov+1):(2*ncov)] + x[2,2,,1:ncov] + u


for (j in (ncov+1):(2*ncov)) {
  xobs[,j] <- x[cbind(z[,1]+1,1,1:ntot,j)]
}
pz <- 1/(1+exp(-(a0 + 2.0 * z[,1] + xobs[,(ncov+1):(2*ncov)] %*% a)))
z[,2] <- (runif(ntot) < pz)

# 3. intervention:

x[1,1,,(2*ncov+1):(3*ncov)] <- rmvnorm(ntot, rep(0.0, ncov), (xdiffsd^2) * s)
x[1,1,,(2*ncov+1):(3*ncov)] <- x[1,1,,(2*ncov+1):(3*ncov)] + x[1,1,,(ncov+1):(2*ncov)] + u

x[1,2,,(2*ncov+1):(3*ncov)] <- rmvnorm(ntot, rep(coeffc, ncov), (xdiffsd^2) * s)
x[1,2,,(2*ncov+1):(3*ncov)] <- x[1,2,,(2*ncov+1):(3*ncov)] + x[1,2,,(ncov+1):(2*ncov)] + u

x[2,1,,(2*ncov+1):(3*ncov)] <- rmvnorm(ntot, rep(0.0, ncov), (xdiffsd^2) * s)
x[2,1,,(2*ncov+1):(3*ncov)] <- x[2,1,,(2*ncov+1):(3*ncov)] + x[2,1,,(ncov+1):(2*ncov)] + u

x[2,2,,(2*ncov+1):(3*ncov)] <- rmvnorm(ntot, rep(coeffc, ncov), (xdiffsd^2) * s)
x[2,2,,(2*ncov+1):(3*ncov)] <- x[2,2,,(2*ncov+1):(3*ncov)] + x[2,2,,(ncov+1):(2*ncov)] + u

for (j in (2*ncov+1):(3*ncov)) {
  xobs[,j] <- x[cbind(z[,1]+1,z[,2]+1,1:ntot,j)]
}
pz <- 1/(1+exp(-(a0 + 0.0 * z[,1] + 2.0 * z[,2] + xobs[,(2*ncov+1):(3*ncov)] %*% a)))
z[,3] <- (runif(ntot) < pz)
zobs <- rowSums(z)
table(zobs)

# Outcome:

ylong <- NULL
zlong <- NULL
counter <- 1
for (j in 0:1) {        
  for (k in 0:1) {
    for (l in 0:1) {
      py <- 1/(1+exp(-(-1.0 + rowMeans(u) + effect * (j+k+l) + x[j+1,k+1,,(2*ncov+1):(3*ncov)] %*% b)))
      y[j+1,k+1,l+1,] <- (runif(ntot) < py)
      ylong <- c(ylong, y[cbind(j+1,k+1,l+1,1:ntot)])
      zlong <- c(zlong, rep(j+k+l, ntot))
    }
  }
}
yobs <- y[cbind(z[,1]+1,z[,2]+1,z[,3]+1,1:ntot)]
table(yobs)

truemodel <- glm(ylong ~ zlong, family=binomial(link=logit))
summary(truemodel)

coef(truemodel)[2]
vcov(truemodel)[2,2]

y <- y[,,,1:nobs]
yobs <- yobs[1:nobs]
z <- z[1:nobs,]
zobs <- zobs[1:nobs]
x <- x[,,1:nobs,]
xobs <- xobs[1:nobs,]

observed.data = data.frame(cbind(y=yobs, xobs, z))
names(observed.data) = c("Y", paste("X", rep(1:3, each = 5), 1:5, sep=""), paste("Z",1:3, sep=""))
write.csv2(observed.data, file = "saerela.csv")
