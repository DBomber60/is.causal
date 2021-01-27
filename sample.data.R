# simulate data like Saaerela
# msmsim('./results0.txt', from=1, to=2, nobs=2000, ncov=5, coeffa=0.3, coeffb=0.3, coeffc=-0.75, effect=-0.25, a0=-1.0, dobayes=TRUE, doboot=TRUE)
library(mvtnorm)
source('mcmc.samples.g.R')
set.seed(1)

# start with a single iteration of simulating data, generate bayesian g formula results
# msmsim = function(outfile, from=1, to=1000, nobs=500, ncov=5, coeffa, coeffb, coeffc, effect, a0, nboot=100) {
# 2. generate parametric g-formula

ncov = 5 # number of (time-varying) covariates
a0=-1
coeffa=0.3; coeffb=0.3; coeffc=-0.75; effect=-0.25
nobs=500
ntot = 1000 # total number of 'true' observations generated (observed is a subset of this)

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

# true 'causal' parameter

predict(truemodel, newdata = data.frame(zlong = 3), type = "response") - 
  predict(truemodel, newdata = data.frame(zlong = 0), type = "response")

y <- y[,,,1:nobs]
yobs <- yobs[1:nobs]
z <- z[1:nobs,]
zobs <- zobs[1:nobs]
x <- x[,,1:nobs,]
xobs <- xobs[1:nobs,]

observed.data = data.frame(cbind(y=yobs, xobs, z))
names(observed.data) = c("Y", paste("X", rep(1:3, each = 5), 1:5, sep=""), paste("Z",1:3, sep=""))
#write.csv(observed.data, file = "saerela.csv")

# bayesian.gf(dat = observed.data, nIter = 100)

# IPT weights:

tmodel1 <- glm(z[,1] ~ xobs[,1:ncov], family=binomial(link=logit))
tmodel2 <- glm(z[,2] ~ xobs[,(ncov+1):(2*ncov)] + z[,1], family=binomial(link=logit))
tmodel3 <- glm(z[,3] ~ xobs[,(2*ncov+1):(3*ncov)] + z[,1] + z[,2], family=binomial(link=logit))
lp1 <- cbind(1.0, xobs[,1:ncov]) %*% coef(tmodel1)
lp2 <- cbind(1.0, xobs[,(ncov+1):(2*ncov)], z[,1]) %*% coef(tmodel2)
lp3 <- cbind(1.0, xobs[,(2*ncov+1):(3*ncov)], z[,1], z[,2]) %*% coef(tmodel3)
pt <- (exp(z[,1] * lp1)/(1+exp(lp1))) * 
  (exp(z[,2] * lp2)/(1+exp(lp2))) *
  (exp(z[,3] * lp3)/(1+exp(lp3)))

smodel1 <- glm(z[,1] ~ 1, family=binomial(link=logit))
smodel2 <- glm(z[,2] ~ z[,1], family=binomial(link=logit))
smodel3 <- glm(z[,3] ~ z[,1] + z[,2], family=binomial(link=logit))
lp1 <- cbind(rep(1.0, nobs)) %*% coef(smodel1)
lp2 <- cbind(1.0, z[,1]) %*% coef(smodel2)
lp3 <- cbind(1.0, z[,1], z[,2]) %*% coef(smodel3)
sc <- (exp(z[,1] * lp1)/(1+exp(lp1))) * 
  (exp(z[,2] * lp2)/(1+exp(lp2))) *
  (exp(z[,3] * lp3)/(1+exp(lp3)))      

iptw <- 1.0/pt
iptws <- sc * iptw

summary(glm(yobs ~ rowSums(z), family = quasibinomial(), weights = iptws))



