#Multiple imputation
#Wide format
#Integer, binary, categorical variables
#Multilevel

library(tidyverse)
library(mice)
dat <- read_csv("C:/Users/akihi/Downloads/Sample_data.csv",
                locale = locale(encoding = "SHIFT-JIS"),
                col_types = cols(
                  id = col_double(),
                  age = col_double(),
                  gender = col_factor(),
                  wbc = col_double(),
                  eosi_p = col_double(),
                  sbp = col_double(),
                  dbp = col_double(),
                  bun = col_double(),
                  rr = col_double(),
                  ams = col_factor(),
                  hr = col_double(),
                  death = col_factor(),
                  adl = col_factor(),
                  hospital = col_factor()
                ),
                guess_max = 1500, #default: 1000
                na = "NA") #Look at the errors
dat %>% glimpse()
md.pattern(dat)
###################################################################################################################################
#1. Generate m complete data sets.(mice -> mids {a multiply imputed data set})
#2. Analyze separately using any standard complete-data technique. (with -> mira {a multiple imputed repeated analysis})
#3. Integrate into an overall set of results using combination rules. (pool -> mipo {a multiply imputed pooled outcome} )
###################################################################################################################################
###First step: creating mids object###
#mice: default "pmm", "logreg (Bayesian logistic regression)", "polyreg (Bayesian polytomous regression)"
#maxit: iterations of the Gibbs sampler
dat1 <- mice(dat, maxit = 0)
dat1$method
dat1$predictorMatrix
dat2 <- mice(dat, m = 10, maxit = 20, printFlag = FALSE, seed = 1234) #maxit 20 is enough.
#Spaghetti plot
#Check intermingled spaghetti
plot(dat2)
#Then m=100 imputation
dat_i <- mice(dat, m = 100, maxit = 20, print = FALSE, seed = 1234)
plot(dat_i)
#Standard imputation vs multilevel imputation
#Only support two-level
###Second step###
com <- complete(dat_i, 1)
fit1 <- lm(wbc ~ age, com)
s <- summary(fit)
s$coefficients
#Manual approach
est <- se <- vector(length = dat_i$m, mode = "list")
for (i in 1:dat_i$m) {
  com <- complete(dat_i, i)
  fit2 <- lm(wbc ~ age, com)
  s <- summary(fit2)
  est[[i]] <- s$coefficients[2, 1]
  se[[i]] <- s$coefficients[2, 2]
  }
library(mitools)
miinf <- MIcombine(est, se)
print(miinf)
miinf$coefficients + c(-1, 1) * qt(0.975, miinf$df) * miinf$variance
#Using with & pool
#Dy default, correct for small samples
mira <- with(dat_i, lm(wbc ~ age))
result <- summary(pool(mira))
print(result[, 1:5])








