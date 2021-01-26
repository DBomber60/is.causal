dat = read.csv("saerela.csv")
library(boot)

# function that takes data and indices as input and generates an estimate of the mean outcome under 
# always treat - mean outcome under never treat strategy

standard = function(data, indices) {
  d = data[indices,]
  #d$interv = -1
  # IPT weights:
  
  tmodel1 <- glm(Z1 ~ X11 + X12 + X13 + X14 + X15, data = d, family=binomial(link=logit))
  tmodel2 <- glm(Z2 ~ Z1 + X21 + X22 + X23 + X24 + X25, data = d, family=binomial(link=logit))
  tmodel3 <- glm(Z3 ~ Z2 + X31 + X32 + X33 + X34 + X35, data = d, family=binomial(link=logit))
  p1 <- dbinom(d$Z1, 1, prob = (predict(tmodel1, type = "response")))
  p2 <- dbinom(d$Z2, 1, prob = (predict(tmodel2, type = "response")))
  p3 <- dbinom(d$Z3, 1, prob = (predict(tmodel3, type = "response")))
  pt <- p1 * p2 * p3
  
  smodel1 <- glm(Z1 ~ 1, data = d, family=binomial(link=logit))
  smodel2 <- glm(Z2 ~ Z1, data = d, family=binomial(link=logit))
  smodel3 <- glm(Z3 ~ Z1 + Z2, data = d, family=binomial(link=logit))
  sp1 <- dbinom(d$Z1, 1, prob = (predict(smodel1, type = "response")))
  sp2 <- dbinom(d$Z2, 1, prob = (predict(smodel2, type = "response")))
  sp3 <- dbinom(d$Z3, 1, prob = (predict(smodel3, type = "response")))
  sc <- sp1 * sp2 * sp3
  
  # never treat

  # always treat
  iptw <- 1.0/pt
  iptws <- sc * iptw
  d$s = with(d, Z1 + Z2 + Z3)
  
  fit = (glm(Y ~ s, data = d, weights = iptws, family = quasibinomial))
  
  ey111 = predict(fit, newdata = data.frame(s=3), type = "response")
  ey000 = predict(fit, newdata = data.frame(s=0), type = "response")
  
  return(ey111-ey000)
}

#results <- boot(data = dat, statistic = standard, R = 100)

#mean(results$t)

