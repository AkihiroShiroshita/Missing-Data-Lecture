###Missing data lecture with R
###Data preparation
library(tidyverse)
library(norm)
url <- "https://cran.r-project.org/src/contrib/Archive/norm2/norm2_2.0.3.tar.gz"
pkgFile <- "norm2_2.0.3.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(lavaan)
library(norm2)
library(jomo)
library(mitools)
library(mice)
dat <- read_csv("C:/Users/akihi/Downloads/MissingDataLecture/Sample_data.csv",
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
                  ams = col_logical(),
                  hr = col_double(),
                  death = col_factor(),
                  adl = col_factor(),
                  hospital = col_factor()
                ),
                guess_max = 1500, #default: 1000
                na = "NA") #Look at the errors
dat %>% glimpse()
md.pattern(dat)
###Missing at random
###EM algorithm
dat_matrix <- data.matrix(dat)
s <- prelim.norm(dat_matrix) #Convert to a matrix
emout <- em.norm(s) #Pre-process
output <- getparam.norm(s, emout) #EM algorithm
sigma <- output[[2]] #Covariance matrix
colnames(sigma) <- c("id", "age", "gender", "wbc", "eosi_p", "sbp", "dbp", "bun", "ams", "hr", "death", "adl", "hospital")
sigma
cov2cor(sigma)
means <- output[[1]]
names(means) <- c("id", "age", "gender", "wbc", "eosi_p", "sbp", "dbp", "bun", "ams", "hr", "death", "adl", "hospital")
means
model <- "death ~ 1+age+gender+eosi_p+sbp+bun+ams+hr+adl"
fit <- sem(model, sample.cov = sigma,
           sample.nobs = 1237,
           sample.mean = means)
summary(fit)
#Then calculate SE with Bootstrap
#But this method is wrong.
###FIML
dat <- dat %>%
  mutate_all(.funs = ~as.numeric(.))
model_f <- "death ~ age+gender+eosi_p+sbp+bun+ams+hr+adl"
fit_f <- sem(model_f, data = dat, missing = "fiml")
summary(fit_f)
#This method is also wrong.
###MCMC
outjomo <- jomo1(dat)
outjomo <- subset(outjomo, Imputation>0)
mi_list <- imputationList(split(outjomo, outjomo$Imputation))
mi_results <- with(mi_list, lm(death ~ age+gender+eosi_p+bun+ams+hr+adl))
summary(pool(as.mira(mi_results)))

