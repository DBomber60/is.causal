dat = read.csv("saerela.csv")


# IPT weights:

tmodel1 <- glm(Z1 ~ X11 + X12 + X13 + X14 + X15, data = dat, family=binomial(link=logit))
tmodel2 <- glm(Z2 ~ Z1 + X21 + X22 + X23 + X24 + X25,data = dat, family=binomial(link=logit))
tmodel3 <- glm(Z3 ~ Z2 + X31 + X32 + X33 + X34 + X35, data = dat, family=binomial(link=logit))
p1 <- dbinom(dat$Z1, 1, prob = (predict(tmodel1, type = "response")))
p2 <- dbinom(dat$Z2, 1, prob = (predict(tmodel2, type = "response")))
p3 <- dbinom(dat$Z3, 1, prob = (predict(tmodel3, type = "response")))
pt <- p1 * p2 * p3

smodel1 <- glm(Z1 ~ 1, data = dat, family=binomial(link=logit))
smodel2 <- glm(Z2 ~ Z1, data = dat, family=binomial(link=logit))
smodel3 <- glm(Z3 ~ Z1 + Z2, data = dat, family=binomial(link=logit))
sp1 <- dbinom(dat$Z1, 1, prob = (predict(smodel1, type = "response")))
sp2 <- dbinom(dat$Z2, 1, prob = (predict(smodel2, type = "response")))
sp3 <- dbinom(dat$Z3, 1, prob = (predict(smodel3, type = "response")))
sc <- sp1 * sp2 * sp3

iptw <- 1.0/pt
iptws <- sc * iptw
dat$s = with(dat, Z1 + Z2 + Z3)

summary(glm(Y ~ s, data = dat, weights = iptws, family = quasibinomial))

