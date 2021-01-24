# sample from posteriors
library(rstan)
library(tidyverse)

obs = read.csv("saerela.csv")
n = nrow(obs)

# gamma parameters 
stg.z1.n = stan_glm(Z1 ~ X11 + X12 + X13 + X14 + X15, data = obs, family = binomial)
stg.z2.n = stan_glm(Z2 ~ Z1 + X21 + X22 + X23 + X24 + X25, data = obs, family = binomial)
stg.z3.n = stan_glm(Z3 ~ Z2 + X31 + X32 + X33 + X34 + X35, data = obs, family = binomial)

# alpha parameters
stg.z1.d = stan_glm(Z1 ~ 1, data = obs, family = binomial)
stg.z2.d = stan_glm(Z2 ~ Z1, data = obs, family = binomial)
stg.z3.d = stan_glm(Z3 ~ Z2 + Z1, data = obs, family = binomial)

# numerator
z1mod.sims = as.matrix(stg.z1.n) # 4k x 1
z2mod.sims = as.matrix(stg.z2.n) # 4k x 2
z3mod.sims = as.matrix(stg.z3.n) # 4k x 2

# 
z1mod0 = as.matrix(stg.z1.d) # 4k x 1
z2mod0 = as.matrix(stg.z2.d) # 4k x 2
z3mod0 = as.matrix(stg.z3.d) # 4k x 2


z1probs.n = plogis(outer(rep(1, n), array(z1mod0)))
z1p = apply(z1probs.n, 2, function(x) dbinom(obs$Z1, 1, prob = x))

z2probs = plogis( as.matrix(cbind(1, select(obs, Z1))) %*% t(z2mod0) ) 
z2p = apply(z2probs, 2, function(x) dbinom(obs$Z2, 1, prob = x))

z3probs = plogis( as.matrix(cbind(1, select(obs, Z2, Z1))) %*% t(z3mod0) ) 
z3p = apply(z3probs, 2, function(x) dbinom(obs$Z3, 1, prob = x))

num = rowMeans(z1p * z2p * z3p) # marginal treatment model

z1probs = plogis( as.matrix(cbind(1, select(obs, X11, X12, X13, X14, X15))) %*% t(z1mod.sims))
z1d = apply(z1probs, 2, function(x) dbinom(obs$Z1, 1, prob = x))

z2probs = plogis( as.matrix(cbind(1, select(obs, Z1, X21, X22, X23, X24, X25))) %*% t(z2mod.sims) )
z2d = apply(z2probs, 2, function(x) dbinom(obs$Z2, 1, prob = x))

z3probs = plogis( as.matrix(cbind(1, select(obs, Z2, X31, X32, X33, X34, X35))) %*% t(z3mod.sims) )
z3d = apply(z3probs, 2, function(x) dbinom(obs$Z3, 1, prob = x))

prod = z1d * z2d * z3d
den = rowMeans(prod) # conditional treatment model

wbayes = num/den
hist(wbayes, breaks = 100, xlim = c(0,10))


smalldat = select(obs, Y, Z1, Z2, Z3) 
smalldat$s = with(smalldat, Z1 + Z2 + Z3)

z = rgamma(n, shape = 1)
q = z/sum(z)


summary(glm(Y~s, data = smalldat, weights = wbayes, family = quasibinomial))
