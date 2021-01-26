# parametric g-formula
library(gfoRmula)
library(tidyverse)

datg = read.csv("saerela.csv")[,-1]
head(datg)

names(datg)[2:16] = paste(letters[1:5], rep(1:3, each = 5), sep="")
#datg$id = 1:nrow(datg)

formatted = datg %>% pivot_longer(c(!Y), names_to = c(".value", "set"),
                     names_pattern = "(.)(.)")

formatted$id = rep(1:nrow(datg), each = 3) # to-do.. this automatically (using regex)
formatted$set = as.numeric(formatted$set) - 1
head(formatted)


# use gfoRmula package
id = 'id'
time_name = 'set'
covnames = c(letters[1:5], 'Z')
outcome_name = 'Y'
covtypes = c(rep('normal',5), 'binary')
histories = c(lagged)
histvars = list(c(letters[1:5], 'Z'))
covparams = list(covmodels=c(a ~ lag1_a + lag1_b+ lag1_c + lag1_d + lag1_e + lag1_Z,
                             b ~ lag1_a + lag1_b+ lag1_c + lag1_d + lag1_e + lag1_Z,
                             c ~ lag1_a + lag1_b+ lag1_c + lag1_d + lag1_e + lag1_Z,
                             d ~ lag1_a + lag1_b+ lag1_c + lag1_d + lag1_e + lag1_Z,
                             e ~ lag1_a + lag1_b+ lag1_c + lag1_d + lag1_e + lag1_Z,
                             Z ~ a + b+ c + d + e + lag1_Z))
ymodel = Y ~ Z + lag1_Z + lag2_Z + a + b + c + d + e
intvars = list('Z', 'Z')
interventions = list( list(c(static, rep(0,3))),
                      list(c(static, rep(1,3))) )  
int_descript = c('Never', 'Always')

b = gformula_binary_eof(formatted, id=id, time_name = time_name, covnames = covnames, covtypes = covtypes,
                   covparams = covparams, histvars = histvars, histories = histories, outcome_name = outcome_name,
                   ymodel = ymodel, intvars = intvars, interventions = interventions, int_descript = int_descript,
                   seed = 1, nsamples = 10, ref_int = 2)

# 0.30787258/ 0.03070607
lower = b$result$`MD lower 95% CI`[2] # lower
upper = b$result$`MD upper 95% CI`[2]
mean = b$result$`Mean difference`[2]


