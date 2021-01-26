# bayesian g formula
library(tidyverse)
library(rstanarm)

dat = read.csv("saerela.csv")[,-1]

# covariate 1 model
# get data in the right format to do one of these


names(dat)[2:16] = paste(letters[1:5], rep(1:3, each = 5), sep="")


formatted = datg %>% pivot_longer(c(!Y), names_to = c(".value", "set"),
                                  names_pattern = "(.)(.)")
head(formatted)

formatted$id = rep(1:500, each = 3) # figure out how to do this automatically

f2 = formatted %>% group_by(id) %>% mutate(alag = lag(a), 
                                           blag = lag(b),
                                           clag = lag(c),
                                           dlag = lag(d),
                                           elag = lag(e),
                                           zlag = lag(Z)) %>% ungroup() %>% filter(set > 1)

#amod = stan_glmer(a ~ (1|id) + alag + zlag, data = f2, family = gaussian())
#bmod = stan_glmer(b ~ (1|id) + blag + zlag, data = f2, family = gaussian())
#cmod = stan_glmer(c ~ (1|id) + clag + zlag, data = f2, family = gaussian())
#dmod = stan_glmer(d ~ (1|id) + dlag + zlag, data = f2, family = gaussian())
#emod = stan_glmer(e ~ (1|id) + elag + zlag, data = f2, family = gaussian())

# outcome model
dat$treat = with(dat, Z1 + Z2 + Z3)
ymod = stan_glm(Y ~ a3 + b3 + c3 + d3 + e3 + treat, data = dat, family = binomial())

nIter = 100

t.model = c(1, 0)

outcomes = array(NA, dim = c(nIter,2))

for(j in 1:length(t.model)) {
  for (it in 1:nIter) {
    print(it)
  
    pmod = 3
    a2.treat = (rnorm(n, mean = as.matrix(cbind(1, select(dat, a1), Z1 = t.model[j])) %*% as.matrix(amod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(amod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    b2.treat = (rnorm(n, mean = as.matrix(cbind(1, select(dat, b1), Z1 = t.model[j])) %*% as.matrix(bmod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(bmod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    c2.treat = (rnorm(n, mean = as.matrix(cbind(1, select(dat, c1), Z1 = t.model[j])) %*% as.matrix(cmod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(cmod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    d2.treat = (rnorm(n, mean = as.matrix(cbind(1, select(dat, d1), Z1 = t.model[j])) %*% as.matrix(dmod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(dmod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    e2.treat = (rnorm(n, mean = as.matrix(cbind(1, select(dat, e1), Z1 = t.model[j])) %*% as.matrix(emod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(emod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    
    a3.treat = (rnorm(n, mean = as.matrix(cbind(1, a2.treat, Z2 = t.model[j])) %*% as.matrix(amod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(amod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    b3.treat = (rnorm(n, mean = as.matrix(cbind(1, b2.treat, Z2 = t.model[j])) %*% as.matrix(bmod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(bmod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    c3.treat = (rnorm(n, mean = as.matrix(cbind(1, c2.treat, Z2 = t.model[j])) %*% as.matrix(cmod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(cmod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    d3.treat = (rnorm(n, mean = as.matrix(cbind(1, d2.treat, Z2 = t.model[j])) %*% as.matrix(dmod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(dmod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    e3.treat = (rnorm(n, mean = as.matrix(cbind(1, e2.treat, Z2 = t.model[j])) %*% as.matrix(emod)[it, 1:3], sd = 0.25)) +
      rnorm(n, mean = as.matrix(emod)[it, ((pmod+1):(pmod+n))], sd = 0.1 )
    
    
    y.outcome = as.matrix(cbind(1, a3.treat, b3.treat, c3.treat, d3.treat, e3.treat, treat=3*t.model[j])) %*% as.matrix(ymod)[it,]
    y.bin = rbinom(n, 1, prob = plogis(y.outcome))
    outcomes[it,j] = sum(y.bin)/n

  }
}

hist(outcomes[,2] - outcomes[,1])





