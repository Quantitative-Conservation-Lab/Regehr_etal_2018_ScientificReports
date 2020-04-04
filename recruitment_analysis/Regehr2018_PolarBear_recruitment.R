## 
## 
## 
## Integrated Population Modeling Provides the First Empirical Estimates of Vital Rates
## and Abundance for Polar Bears in the Chukchi Sea
##
## E.V. Regehr, N.J. Hostetter, R.R. Wilson, K.D. Rode, M. St. Martin, and S.J. Converse
##
## Scientific Reports
## October 2018
## 
##
## SUPPLEMENT: DATA AND JAGS MODELS FOR RECRUITMENT ANALYSIS (YEARLINGS [C1S] PER ADULT FEMALE)
##             PERFORMED SEPARATELY FROM OTHER ANALYSES
## 
## REQUIRED FILES:
##   1. Regehr2018_C1perAF_data.txt                  # Data
##   2. pois_glm_XXXX.txt 			                     # JAGS models (provided in this script), where XXXX is a descriptor of each model
##
## DISCLAIMER: The authors provide no guarantee regarding the completeness or functionality of these data and programs, 
## and are not responsible for any consequences of their use. 

#clear all objects in workspace
rm(list=ls())

#required packages
library(jagsUI)

#load data
dat <- read.table( "Regehr2018_C1perAF_data.txt", header = TRUE, sep = " " )

#format data
year<-rownames(table(dat$year))
nAF<-as.vector(table(dat$year))   #annual number of AF's
nC1<-NA                           #annual number of c1 cubs (placeholder)
 for(t in 1:length(year)){ 
  nC1[t]<-sum(dat$nC1[which(dat$year==year[t])])
 }
#standardize year
yearS<-as.numeric(year)-min(as.numeric(year))+1


##BUGS model for linear trend
sink("pois_glm_linear.txt")
cat("
model {
 #priors
 beta0~dnorm(0,0.01)
 beta1~dnorm(0,0.01)

 #likelihood
 for(t in 1:nyears){
  nC1[t]~dpois(lambda[t]*nAF[t])
  lambda[t]<- exp(beta0 + beta1*year[t])
 }
}#model
",fill = TRUE)
sink()

#data
jags.data <- list(nC1=nC1, nAF=nAF, year=yearS, nyears=length(yearS))

#Initial values
inits <- function() list(beta0=runif(1,-1,1), beta1=runif(1, -1, 1))

#parameters monitored
params <- c("lambda", "beta0", "beta1")

#MCMC settings
#nc=3; nAdapt=1000; nb=2000; ni=5000+nb; nt=1
nc=3; nAdapt=10000; nb=20000; ni=1500000+nb; nt=100 #IPM settings

#call JAGS from R
out <- jags(jags.data, inits, params, "pois_glm_linear.txt",
  n.adapt=nAdapt,  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
  parallel=T)

print(out,2)


##BUGS model for time constant instead of trend
sink("pois_glm_constant.txt")
cat("
model {
  beta0~dnorm(0,0.01)
  lambda<- exp(beta0)

 for(t in 1:nyears){
  nC1[t]~dpois(lambda*nAF[t])
 }
}#model
",fill = TRUE)
sink()

#initial values
inits <- function() list(beta0=runif(1,-1,1))

#parameters monitored
params <- c("lambda", "beta0")

#MCMC settings
#nc=3; nAdapt=1000; nb=2000; ni=5000+nb; nt=1
nc=3; nAdapt=10000; nb=20000; ni=1500000+nb; nt=100 #IPM settings

#call JAGS from R
outB <- jags(jags.data, inits, params, "pois_glm_constant.txt",
  n.adapt=nAdapt,  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
  parallel=T)

print(outB)


##
##
##
## END 
##
##
##