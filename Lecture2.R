#Propensity score analysis
#Complete case
#Remove "ID"
#"hospitalterm" is recognized as "normally distributed" just for simplicity
library(tidyverse)
library(survey)
library(mitools)
library(survey)
library(lattice)
library(twang)
library(mice)
dat <- read_csv("C:/Users/akihi/Downloads/MissingDataLecture/Sample_data2.csv",
                locale = locale(encoding = "SHIFT-JIS"),
                col_types = cols(
                  id = col_double(),
                  age = col_double(),
                  steroid = col_factor(),
                  gender = col_factor(),
                  sbp = col_double(),
                  dbp = col_double(),
                  bun = col_double(),
                  ams = col_factor(),
                  hr = col_double(),
                  death = col_factor(),
                  adl = col_factor(),
                  hospital = col_factor(),
                  care = col_factor(),
                  hospitalterm = col_double()
                ),
                guess_max = 1500, #default: 1000
                na = "NA")
dat %>% glimpse()
#dat <- dat %>%
#  na.omit()
long.imputation <- c()
predictor.selection <- quickpred(dat,
                                 exclude=c("id"))
imputation <- mice(dat,
                   m=5,
                   predictorMatrix = predictor.selection)
long.imputation <- rbind(long.imputation,complete(imputation, action="long"))
#####################################################
#1. Propensity score estimation                 #####
#2. Propensity score method implementation      #####
#3. Covariate balance evaluation                #####
#4. Treatment effect estimation                 #####
#6. Sensitivity analysis                        #####
#####################################################
#1. Propensity score estimation
imputation1 = subset(long.imputation, subset=.imp==1)
allimputations <- imputationList(list(
  subset(long.imputation, subset=.imp==1),
  subset(long.imputation, subset=.imp==2),
  subset(long.imputation, subset=.imp==3),
  subset(long.imputation, subset=.imp==4),
  subset(long.imputation, subset=.imp==5)))
Names <- colnames(dat)
psFormula <- paste(Names, collapse="+")
psFormula <- formula(paste("steroid~",psFormula, sep=""))
print(psFormula)
surveyDesign1 <- svydesign(ids=~id,
                           strata=~hospital,
                           data = imputation1,
                           nest=T)
surveyDesignAll <- svydesign(ids=~id,
                             strata=~hospital,
                             data = allimputations,
                             nest=T)
psModel1 <- svyglm(psFormula,
                   design=surveyDesign1,
                   family=binomial)
pScores <- fitted(psModel1)
imputation1$pScores <- pScores
psModelAll <- with(surveyDesignAll, svyglm(psFormula, family=binomial))
pScoresAll <- sapply(psModelAll, fitted)
pScores <- apply(pScoresAll,1,mean)
allimputations <- update(allimputations, pScores = pScores)
#common support
comp <- complete(imputation, action="long")
bwplot(pScores ~ steroid,
       data = comp,
       ylab = "Propensity Scores",
       xlab = "Treatment",
       auto.key = TRUE)
#2. Propensity score method implementation
#PS weight
comp$weightATT <- with(comp,
                       ifelse(steroid==1, 1, pScores/(1-pScores)))
with(comp, by(weightATT,steroid,summary))
comp$weightATE <- with(comp,
                       ifelse(steroid==1, 1/pScores, 1/(1-pScores)))
with(comp, by(weightATE,steroid,summary))
#Truncation
comp$weightATETruncated <- with(comp,
                                ifelse(weightATE > quantile(weightATE, 0.99),
                                                   quantile(weightATE, 0.99),weightATE))
with(comp, by(weightATETruncated,steroid,summary))
#Stabilized weight
comp$C <- with(comp,ifelse(steroid==1,pScores,1-pScores))
surveyDesign <- svydesign(ids=~id,
                          strata=~hospital,
                          data = comp,
                          nest=T)
constants <- svyby(~C, by=~steroid, design=surveyDesign, FUN=svymean)
comp$stabilizedWeightATE <- ifelse(comp$steroid==1,
                                   constants[1,2]/comp$C,
                                   constants[2,2]/comp$C)
with(comp, by(stabilizedWeightATE,steroid,summary))
#Balance check
balanceTable <- bal.stat(comp, vars= Names,
                         treat.var = "steroid",
                         w.all = comp$weightATETruncated, get.ks=F,
                         sampw = 1,
                         estimand="ATE", multinom=F)
balanceTable <- balanceTable$results
round(balanceTable,3)
#Estimating treatment effects
surveyDesignLS <- svydesign(ids=~id,
                            strata=~hospital,
                            weights=~weightATETruncated,
                            data = comp,
                            nest=T)
surveyDesignLSBoot <- as.svrepdesign(surveyDesignLS, type=c("bootstrap"), replicates=1000)
weightedMeans <- svyby(formula=~hospitalterm,
                        by=~steroid,
                        surveyDesignLSBoot,
                        FUN=svymean,
                        covmat=TRUE)
weightedMeans
weightedVars <- svyby(formula=~hospitalterm,
                      by=~steroid,
                      design=surveyDesignLSBoot,
                      FUN=svyvar,
                      covmat=TRUE)
weightedVars
Effect <- svycontrast(weightedMeans, contrasts=c(-1,1))
Effect
#Weighted least squares
outcomeModel <- svyglm(hospitalterm ~ steroid, surveyDesignLS)
summary(outcomeModel)
#However, is the method described so far truly correct?
allimputations <- update(allimputations,
                         weightATT = ifelse(steroid==1, 1, pScores/(1-pScores)))
surveyDesignMI <- svydesign(ids=~id,
                            strata=~hospital,
                            weights=~weightATT,
                            data = allimputations,
                            nest=T)
outcomeModelMI <- with(surveyDesignMI, svyglm(hospitalterm~steroid))
resultsModelMI <- MIcombine(outcomeModelMI)
summary(resultsModelMI)
#Doubly robust estimation
#Add PS to the data
surveyDesignT <- subset(surveyDesignLSBoot, steroid==1)
surveyDesignC <- subset(surveyDesignLSBoot, steroid==0)
modelT <- svyglm(hospitalterm ~ pScores+I(pScores^2)+I(pScores^3), surveyDesignT)
modelU <- svyglm(hospitalterm ~ pScores+I(pScores^2)+I(pScores^3), surveyDesignC)
Yt1 <- predict(modelT,
               newdata=data.frame(pScores=with(comp,pScores[steroid==1])),
               vcov=TRUE, type="response")
Yt0 <- predict(modelC,newdata=data.frame(pScores=with(comp,pScores[steroid==1])),
               vcov=TRUE, type="response")
diff <- Yt1 - Yt0
ATTdr <- svycontrast(diff, contrasts=rep(1/length(diff),length(diff)))
#Stratification














