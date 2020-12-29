#Multiple imputation
#Wide format
#Integer, binary, categorical variables
#Multilevel
#Set up libraries
library(tidyverse)
library(mice)
url <- "https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz"
pkgFile <- "norm2_2.0.3.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(norm2)
#Read the sample file
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
############################################################################################################################
#1. Generate m complete data sets.(mice -> mids {a multiply imputed data set})                                            ##
#2. Analyze separately using any standard complete-data technique. (with -> mira {a multiple imputed repeated analysis})  ##
#3. Integrate into an overall set of results using combination rules. (pool -> mipo {a multiply imputed pooled outcome} ) ##
############################################################################################################################
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
###Second and third step###
com <- complete(dat_i, 1)
fit1 <- glm(death ~ age, com, family = binomial(link = "logit"))
s <- summary(fit1)
exp(s$coefficients)
#Manual approach
est <- se <- vector(length = dat_i$m, mode = "list")
for (i in 1:dat_i$m) {
  com <- complete(dat_i, i)
  fit2 <- glm(death ~ age, com, family = binomial(link = "logit"))
  s <- summary(fit2)
  est[[i]] <- s$coefficients[2, 1]
  se[[i]] <- s$coefficients[2, 2]
  }
miinf <- miInference(est, se)
print(miinf)
exp(miinf$est)
exp(miinf$std.err)
#Using with & pool
#Dy default, correct for small samples
mira <- with(dat_i, glm(death ~ age,  family = binomial))
result <- summary(pool(mira))
print(result[, 1:5])
##Model diagnostics
imp.wbc <- data.frame(dat_i$imp$wbc)
imp.wbc.long <- reshape(imp.wbc,
                        varying = list(c(paste0("X", 1:100))),
                        direction = "long",
                        v.names = "wbc",
                        times = 1:100,
                        timevar = "IMP",
                        idvar = "id")
sub <- subset(imp.wbc.long, imp.wbc.long$IMP <= 20)
boxplot(wbc  ~  IMP, sub, xlab = "Imputation number", ylab = "WBC", pch = 20, col = "gray80")
hist(sub$wbc, col = gray(0.1, alpha = 0.5), freq = FALSE, xlab = "Distribution of WBC", main = "")
hist(imp.wbc.long$wbc, col = gray(0.8, alpha = 0.5), freq = FALSE, add = TRUE)
legend("topright", legend = c("observed", "imputed"), fill = c(gray(0.1, alpha = 0.5), gray(0.8, alpha = 0.5)))

