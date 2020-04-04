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
## SUPPLEMENT: DATA AND JAGS MODEL FOR THE INTEGRATED POPULATION MODEL 
## 
## REQUIRED FILES:
##   1. Regehr2018_ jags_data.txt                    # Data
##   2. Regehr2018_ state_inits.txt                  # initial values
##   3. Regehr2018_PolarBear_IPM_JAGS_model_GOF.txt  # JAGS model (provided at end of this script).
##
## DISCLAIMER: The authors provide no guarantee regarding the completeness or functionality of these data and programs, 
## and are not responsible for any consequences of their use. 

#clear all objects in workspace
rm(list=ls())

#required packages
library(jagsUI)

# Load data
jags.data <- dget("jags_data.txt")

#Options for first year stage distribution
ssd<-c(0.05,0.04,0.15,0.13,0.09,0.07,0.05,0.04,0.03,0.35)  # SSD used in manuscript (Regehr et al. 2017)
#ssd<-c(0.04,0.03,0.18,0.12,0.09,0.07,0.04,0.03,0.03,0.36) # Alternative SSD (see Supplementary Methods)
#ssd<-c(0.06,0.05,0.11,0.12,0.10,0.08,0.06,0.05,0.04,0.35) # Alternative SSD (see Supplementary Methods) 


#load intial values for state processes (reasonable starting values prevents JAGS from crashing)
state_inits <- dget("state_inits.txt")

#initial values
inits <- function(){list(
 zF=state_inits$zF.inits, zM=state_inits$zM.inits,	
 phiFAD = runif(1,.8,1), phiFSA = runif(1,.6,.8),
 phiMAD = runif(1,.8,1), phiMSA = runif(1,.6,.8),
 phiCub0 = runif(1,.6,.8), phiCub1 = runif(1,.8,1),
 pDet = runif(1,.20,.30), B1 = runif(1,.6,.8), B2 = runif(1,.15,.25),
 W = runif(1,.3,.4), psiF = c(0.60,0.75,0.05,0.90),
 R0p = c(0.20, 0.42, 0.38), Nstate=state_inits$NstateInits
)}#inits

#monitored parameters
params <- c("phiFAD","phiFSA","phiMAD","phiMSA",
 "phiCub0","phiCub1","phiL0","phiL1",			
 "pDet", "pFail", "B1","B2",  "W","psiF","psiM",
 "R0","R1","C0litterSize","C1litterSize",
 "aAFNCO","aAFNCU","aAFC0O","aAFC0U",
 "aAFC1O","aAFC1U","aAFC2O","aAFC2U","aDEAD",
 "Nc0O","Nc0U","Nc1O","Nc1U",
 "Nstate","Nin",
 "CS_N_fit", "CS_N_fitNew",
 "CS_LS_fit", "CS_LS_fitNew",
 "CS_NW_fit", "CS_NW_fitNew"
 )#params

#MCMC values
nc=3; nAdapt=100; nb=10; ni=100+nb; nt=1
#nc=3; nAdapt=10000; nb=20000; ni=1500000+nb; nt=100	


#Run jagsUI
out <- jags(jags.data, inits, params, "Regehr2018_PolarBear_IPM_JAGS_model_GOF.txt",
  n.adapt=nAdapt,  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)

print(out,2)


#Various summary statistics
mean(out$sims.list$CS_LS_fit<out$sims.list$CS_LS_fitNew)  #GOF for litter size data
mean(out$sims.list$CS_N_fit<out$sims.list$CS_N_fitNew)    #GOF for count data
mean(out$sims.list$CS_NW_fit<out$sims.list$CS_NW_fitNew)  #GOF for weaning data




##
##
##
## JAGS MODEL
##
##
##


cat("
model {

### Component I. CMR FEMALE survival, detection, and state transitions
## Multievent model separated by sex.
 for(i in 1:nFem){
 for(t in (firstFem[i]+1):nyrs){
    #State process (draw zF[t] given zF[t-1])
      zF[i,t] ~ dcat(psF[zF[i,t-1], ])
    #Events (draw event[t] given zF[t]) 
      yF[i,t] ~ dcat(poF[zF[i,t], i, t-1,])
 }#t
}#i


# Female State and Observation matrices
# [state from,individual, time, state to]

psF[1,1] <- 0                                 # F2Y/O to F2Y/O
psF[1,2] <- 0                                 # F2Y/O to F2Y/U
psF[1,3] <- phiFSA*psiF[1]                    # F2Y/O to F3Y/O
psF[1,4] <- phiFSA*(1-psiF[1])                # F2Y/O to F3Y/U
psF[1,5] <- 0                                 # F2Y/O to AFNC/O
psF[1,6] <- 0                                 # F2Y/O to AFNC/U
psF[1,7] <- 0                                 # F2Y/O to AFC0/O
psF[1,8] <- 0                                 # F2Y/O to AFC0/U
psF[1,9] <- 0                                 # F2Y/O to AFC1/O
psF[1,10]<- 0                                 # F2Y/O to AFC1/U
psF[1,11]<- 0                                 # F2Y/O to AFC2/O
psF[1,12]<- 0                                 # F2Y/O to AFC2/U
psF[1,13]<- (1-phiFSA)                        # F2Y/O to Dead
psF[2,1] <- 0                                 # F2Y/U to F2Y/O
psF[2,2] <- 0                                 # F2Y/U to F2Y/U
psF[2,3] <- phiFSA*(1-psiF[2])                # F2Y/U to F3Y/O
psF[2,4] <- phiFSA*psiF[2]                    # F2Y/U to F3Y/U
psF[2,5] <- 0				                  # F2Y/U to AFNC/O
psF[2,6] <- 0                                 # F2Y/U to AFNC/U
psF[2,7] <- 0                                 # F2Y/U to AFC0/O
psF[2,8] <- 0                                 # F2Y/U to AFC0/U
psF[2,9] <- 0                                 # F2Y/U to AFC1/O
psF[2,10]<- 0                                 # F2Y/U to AFC1/U
psF[2,11]<- 0                                 # F2Y/U to AFC2/O
psF[2,12]<- 0                                 # F2Y/U to AFC2/U
psF[2,13]<- (1-phiFSA)                        # F2Y/U to Dead
psF[3,1] <- 0                                 # F3Y/O to F2Y/O
psF[3,2] <- 0                                 # F3Y/O to F2Y/U
psF[3,3] <- 0                                 # F3Y/O to F3Y/O
psF[3,4] <- 0                                 # F3Y/O to F3Y/U
psF[3,5] <- phiFSA*psiF[1]                    # F3Y/O to AFNC/O
psF[3,6] <- phiFSA*(1-psiF[1])                # F3Y/O to AFNC/U
psF[3,7] <- 0                                 # F3Y/O to AFC0/O
psF[3,8] <- 0                                 # F3Y/O to AFC0/U
psF[3,9] <- 0                                 # F3Y/O to AFC1/O
psF[3,10]<- 0                                 # F3Y/O to AFC1/U
psF[3,11]<- 0                                 # F3Y/O to AFC2/O
psF[3,12]<- 0                                 # F3Y/O to AFC2/U
psF[3,13]<- (1-phiFSA)                        # F3Y/O to Dead
psF[4,1] <- 0                                 # F3Y/U to F2Y/O
psF[4,2] <- 0                                 # F3Y/U to F2Y/U
psF[4,3] <- 0                                 # F3Y/U to F3Y/O
psF[4,4] <- 0                                 # F3Y/U to F3Y/U
psF[4,5] <- phiFSA*(1-psiF[2])                # F3Y/U to AFNC/O
psF[4,6] <- phiFSA*psiF[2]                    # F3Y/U to AFNC/U
psF[4,7] <- 0                                 # F3Y/U to AFC0/O
psF[4,8] <- 0                                 # F3Y/U to AFC0/U
psF[4,9] <- 0                                 # F3Y/U to AFC1/O
psF[4,10]<- 0                                 # F3Y/U to AFC1/U
psF[4,11]<- 0                                 # F3Y/U to AFC2/O
psF[4,12]<- 0                                 # F3Y/U to AFC2/U
psF[4,13]<- (1-phiFSA)                        # F3Y/U to Dead
psF[5,1] <- 0                                 # AFNC/O to F2Y/O
psF[5,2] <- 0                                 # AFNC/O to F2Y/U
psF[5,3] <- 0                                 # AFNC/O to F3Y/O
psF[5,4] <- 0                                 # AFNC/O to F3Y/U
psF[5,5] <- phiFAD*psiF[1]*(1-B1)             # AFNC/O to AFNC/O
psF[5,6] <- phiFAD*(1-psiF[1])*(1-B1)         # AFNC/O to AFNC/U
psF[5,7] <- phiFAD*psiF[3]*(B1)               # AFNC/O to AFC0/O
psF[5,8] <- phiFAD*(1-psiF[3])*(B1)           # AFNC/O to AFC0/U
psF[5,9] <- 0                                 # AFNC/O to AFC1/O
psF[5,10]<- 0                                 # AFNC/O to AFC1/U
psF[5,11]<- 0                                 # AFNC/O to AFC2/O
psF[5,12]<- 0                                 # AFNC/O to AFC2/U
psF[5,13]<- (1-phiFAD)                        # AFNC/O to Dead
psF[6,1] <- 0                                 # AFNC/U to F2Y/O
psF[6,2] <- 0                                 # AFNC/U to F2Y/U
psF[6,3] <- 0                                 # AFNC/U to F3Y/O
psF[6,4] <- 0                                 # AFNC/U to F3Y/U
psF[6,5] <- phiFAD*(1-psiF[2])*(1-B1)         # AFNC/U to AFNC/O
psF[6,6] <- phiFAD*psiF[2]*(1-B1)             # AFNC/U to AFNC/U
psF[6,7] <- phiFAD*(1-psiF[4])*(B1)           # AFNC/U to AFC0/O
psF[6,8] <- phiFAD*psiF[4]*(B1)               # AFNC/U to AFC0/U
psF[6,9] <- 0                                 # AFNC/U to AFC1/O
psF[6,10]<- 0                                 # AFNC/U to AFC1/U
psF[6,11]<- 0                                 # AFNC/U to AFC2/O
psF[6,12]<- 0                                 # AFNC/U to AFC2/U
psF[6,13]<- (1-phiFAD)                        # AFNC/U to Dead
psF[7,1] <- 0                                     # AFC0/O to F2Y/O
psF[7,2] <- 0                                     # AFC0/O to F2Y/U
psF[7,3] <- 0                                     # AFC0/O to F3Y/O
psF[7,4] <- 0                                     # AFC0/O to F3Y/U
psF[7,5] <- phiFAD*psiF[1]*(1-B2)*(1-phiL0)       # AFC0/O to AFNC/O
psF[7,6] <- phiFAD*(1-psiF[1])*(1-B2)*(1-phiL0)   # AFC0/O to AFNC/U
psF[7,7] <- phiFAD*psiF[3]*(B2)*(1-phiL0)         # AFC0/O to AFC0/O
psF[7,8] <- phiFAD*(1-psiF[3])*(B2)*(1-phiL0)     # AFC0/O to AFC0/U
psF[7,9] <- phiFAD*psiF[1]*(phiL0)                # AFC0/O to AFC1/O
psF[7,10]<- phiFAD*(1-psiF[1])*(phiL0)            # AFC0/O to AFC1/U
psF[7,11]<- 0                                     # AFC0/O to AFC2/O
psF[7,12]<- 0                                     # AFC0/O to AFC2/U
psF[7,13]<- (1-phiFAD)                            # AFC0/O to Dead
psF[8,1] <- 0                                     # AFC0/U to F2Y/O
psF[8,2] <- 0                                     # AFC0/U to F2Y/U
psF[8,3] <- 0                                     # AFC0/U to F3Y/O
psF[8,4] <- 0                                     # AFC0/U to F3Y/U
psF[8,5] <- phiFAD*(1-psiF[2])*(1-B2)*(1-phiL0)   # AFC0/U to AFNC/O
psF[8,6] <- phiFAD*(psiF[2])*(1-B2)*(1-phiL0)     # AFC0/U to AFNC/U
psF[8,7] <- phiFAD*(1-psiF[4])*(B2)*(1-phiL0)     # AFC0/U to AFC0/O
psF[8,8] <- phiFAD*(psiF[4])*(B2)*(1-phiL0)       # AFC0/U to AFC0/U
psF[8,9] <- phiFAD*(1-psiF[2])*phiL0              # AFC0/U to AFC1/O
psF[8,10]<- phiFAD*(psiF[2])*phiL0                # AFC0/U to AFC1/U
psF[8,11]<- 0                                     # AFC0/U to AFC2/O
psF[8,12]<- 0                                     # AFC0/U to AFC2/U
psF[8,13]<- (1-phiFAD)                            # AFC0/U to Dead
psF[9,1] <- 0                                     # AFC1/O to F2Y/O
psF[9,2] <- 0                                     # AFC1/O to F2Y/U
psF[9,3] <- 0                                     # AFC1/O to F3Y/O
psF[9,4] <- 0                                     # AFC1/O to F3Y/U
psF[9,5] <- phiFAD*psiF[1]*(1-B2)*(1-phiL1) +
              phiFAD*psiF[1]*phiL1*W[1]           # AFC1/O to AFNC/O
psF[9,6] <- phiFAD*(1-psiF[1])*(1-B2)*(1-phiL1) +
              phiFAD*(1-psiF[1])*phiL1*W[1]       # AFC1/O to AFNC/U
psF[9,7] <- phiFAD*psiF[3]*(B2)*(1-phiL1)         # AFC1/O to AFC0/O
psF[9,8] <- phiFAD*(1-psiF[3])*(B2)*(1-phiL1)     # AFC1/O to AFC0/U
psF[9,9] <- 0                                     # AFC1/O to AFC1/O
psF[9,10]<- 0                                     # AFC1/O to AFC1/U
psF[9,11]<- phiFAD*psiF[1]*phiL1*(1-W[1])         # AFC1/O to AFC2/O
psF[9,12]<- phiFAD*(1-psiF[1])*phiL1*(1-W[1])     # AFC1/O to AFC2/U
psF[9,13]<- (1-phiFAD)                            # AFC1/O to Dead
psF[10,1]<- 0                                     # AFC1/U to F2Y/O
psF[10,2]<- 0                                     # AFC1/U to F2Y/U
psF[10,3]<- 0                                     # AFC1/U to F3Y/O
psF[10,4]<- 0                                     # AFC1/U to F3Y/U
psF[10,5]<- phiFAD*(1-psiF[2])*(1-B2)*(1-phiL1)+
              phiFAD*(1-psiF[2])*phiL1*W[1]       # AFC1/U to AFNC/O
psF[10,6]<- phiFAD*psiF[2]*(1-B2)*(1-phiL1)+
              phiFAD*psiF[2]*phiL1*W[1]           # AFC1/U to AFNC/U
psF[10,7] <- phiFAD*(1-psiF[4])*(B2)*(1-phiL1)    # AFC1/U to AFC0/O
psF[10,8] <- phiFAD*psiF[4]*(B2)*(1-phiL1)        # AFC1/U to AFC0/U
psF[10,9] <- 0                                    # AFC1/U to AFC1/O
psF[10,10]<- 0                                    # AFC1/U to AFC1/U
psF[10,11]<- phiFAD*(1-psiF[2])*phiL1*(1-W[1])    # AFC1/U to AFC2/O
psF[10,12]<- phiFAD*(psiF[2])*phiL1*(1-W[1])      # AFC1/U to AFC2/U
psF[10,13]<- (1-phiFAD)                           # AFC1/U to Dead
psF[11,1] <- 0                                    # AFC2/O to F2Y/O
psF[11,2] <- 0                                    # AFC2/O to F2Y/U
psF[11,3] <- 0                                    # AFC2/O to F3Y/O
psF[11,4] <- 0                                    # AFC2/O to F3Y/U
psF[11,5] <- phiFAD*(psiF[1])*(1-B1)              # AFC2/O to AFNC/O
psF[11,6] <- phiFAD*(1-psiF[1])*(1-B1)            # AFC2/O to AFNC/U
psF[11,7] <- phiFAD*(psiF[3])*(B1)                # AFC2/O to AFC0/O
psF[11,8] <- phiFAD*(1-psiF[3])*(B1)              # AFC2/O to AFC0/U
psF[11,9] <- 0                                    # AFC2/O to AFC1/O
psF[11,10]<- 0                                    # AFC2/O to AFC1/U
psF[11,11]<- 0                                    # AFC2/O to AFC2/O
psF[11,12]<- 0                                    # AFC2/O to AFC2/U
psF[11,13]<- (1-phiFAD)                           # AFC2/O to Dead
psF[12,1] <- 0                                    # AFC2/U to F2Y/O
psF[12,2] <- 0                                    # AFC2/U to F2Y/U
psF[12,3] <- 0                                    # AFC2/U to F3Y/O
psF[12,4] <- 0                                    # AFC2/U to F3Y/U
psF[12,5] <- phiFAD*(1-psiF[2])*(1-B1)            # AFC2/U to AFNC/O
psF[12,6] <- phiFAD*(psiF[2])*(1-B1)              # AFC2/U to AFNC/U
psF[12,7] <- phiFAD*(1-psiF[4])*(B1)              # AFC2/U to AFC0/O
psF[12,8] <- phiFAD*(psiF[4])*(B1)                # AFC2/U to AFC0/U
psF[12,9] <- 0                                    # AFC2/U to AFC1/O
psF[12,10]<- 0                                    # AFC2/U to AFC1/U
psF[12,11]<- 0                                    # AFC2/U to AFC2/O
psF[12,12]<- 0                                    # AFC2/U to AFC2/U
psF[12,13]<- (1-phiFAD)                           # AFC2/U to Dead
psF[13,1] <- 0                             # Dead to F2Y/O
psF[13,2] <- 0                             # Dead to F2Y/U
psF[13,3] <- 0                             # Dead to F3Y/O
psF[13,4] <- 0                             # Dead to F3Y/U
psF[13,5] <- 0                             # Dead to AFNC/O
psF[13,6] <- 0                             # Dead to AFNC/U
psF[13,7] <- 0                             # Dead to AFC0/O
psF[13,8] <- 0                             # Dead to AFC0/U
psF[13,9] <- 0                             # Dead to AFC1/O
psF[13,10]<- 0                             # Dead to AFC1/U
psF[13,11]<- 0                             # Dead to AFC2/O
psF[13,12]<- 0                             # Dead to AFC2/U
psF[13,13]<- 1                             # Dead to Dead

# End Female state matrix

# Female observation matrix
# poF[state,ind,time,event type]
for(i in 1:nFem){
for(t in firstFem[i]:(nyrs-1)){
poF[1,i,t,1] <- pDetT[t]               #F2Y/O
poF[1,i,t,2] <- 0
poF[1,i,t,3] <- 0
poF[1,i,t,4] <- 0
poF[1,i,t,5] <- 0
poF[1,i,t,6] <- 0
poF[1,i,t,7] <- 0
poF[1,i,t,8] <- 0
poF[1,i,t,9] <- 0
poF[1,i,t,10]<- 0
poF[1,i,t,11]<- 0
poF[1,i,t,12]<- 0
poF[1,i,t,13]<- 0
poF[1,i,t,14]<- 0
poF[1,i,t,15]<- 0
poF[1,i,t,16]<- 0
poF[1,i,t,17]<- (1-pDetT[t])

poF[2,i,t,1] <- 0                  #F2Y/U
poF[2,i,t,2] <- 0
poF[2,i,t,3] <- 0
poF[2,i,t,4] <- 0
poF[2,i,t,5] <- 0
poF[2,i,t,6] <- 0
poF[2,i,t,7] <- 0
poF[2,i,t,8] <- 0
poF[2,i,t,9] <- 0
poF[2,i,t,10]<- 0
poF[2,i,t,11]<- 0
poF[2,i,t,12]<- 0
poF[2,i,t,13]<- 0
poF[2,i,t,14]<- 0
poF[2,i,t,15]<- 0
poF[2,i,t,16]<- 0
poF[2,i,t,17]<- 1

poF[3,i,t,1] <- 0                  #F3Y/O
poF[3,i,t,2] <- pDetT[t]
poF[3,i,t,3] <- 0
poF[3,i,t,4] <- 0
poF[3,i,t,5] <- 0
poF[3,i,t,6] <- 0
poF[3,i,t,7] <- 0
poF[3,i,t,8] <- 0
poF[3,i,t,9] <- 0
poF[3,i,t,10]<- 0
poF[3,i,t,11]<- 0
poF[3,i,t,12]<- 0
poF[3,i,t,13]<- 0
poF[3,i,t,14]<- 0
poF[3,i,t,15]<- 0
poF[3,i,t,16]<- 0
poF[3,i,t,17]<- (1-pDetT[t])

poF[4,i,t,1] <- 0                  #F3Y/U
poF[4,i,t,2] <- 0
poF[4,i,t,3] <- 0
poF[4,i,t,4] <- 0
poF[4,i,t,5] <- 0
poF[4,i,t,6] <- 0
poF[4,i,t,7] <- 0
poF[4,i,t,8] <- 0
poF[4,i,t,9] <- 0
poF[4,i,t,10]<- 0
poF[4,i,t,11]<- 0
poF[4,i,t,12]<- 0
poF[4,i,t,13]<- 0
poF[4,i,t,14]<- 0
poF[4,i,t,15]<- 0
poF[4,i,t,16]<- 0
poF[4,i,t,17]<- 1

poF[5,i,t,1] <- 0                  #AFNC/O
poF[5,i,t,2] <- 0
poF[5,i,t,3] <- pDetT[t]*(1-Telem[i,t])+pDetT[t]*((Telem[i,t]*pFail[TelemYears[i,t],t]))
poF[5,i,t,4] <- 0
poF[5,i,t,5] <- 0
poF[5,i,t,6] <- 0
poF[5,i,t,7] <- aAFNCO[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[5,i,t,8] <- 0
poF[5,i,t,9] <- 0
poF[5,i,t,10]<- 0
poF[5,i,t,11]<- aAFNCO[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[5,i,t,12]<- aAFNCO[3]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[5,i,t,13]<- 0
poF[5,i,t,14]<- aAFNCO[4]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[5,i,t,15]<- aAFNCO[5]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[5,i,t,16]<- 0
poF[5,i,t,17]<- (1-pDetT[t])*(1-Telem[i,t])+(1-pDetT[t])*((Telem[i,t]*pFail[TelemYears[i,t],t]))

poF[6,i,t,1] <- 0                  #AFNC/U
poF[6,i,t,2] <- 0
poF[6,i,t,3] <- 0
poF[6,i,t,4] <- 0
poF[6,i,t,5] <- 0
poF[6,i,t,6] <- 0        
poF[6,i,t,7] <- 0
poF[6,i,t,8] <- 0
poF[6,i,t,9] <- 0
poF[6,i,t,10]<- 0
poF[6,i,t,11]<- aAFNCU[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])    
poF[6,i,t,12]<- 0
poF[6,i,t,13]<- aAFNCU[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[6,i,t,14]<- aAFNCU[3]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[6,i,t,15]<- 0
poF[6,i,t,16]<- aAFNCU[4]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[6,i,t,17]<- (1-Telem[i,t])+(Telem[i,t]*pFail[TelemYears[i,t],t])         

poF[7,i,t,1] <- 0                  #AFC0/O
poF[7,i,t,2] <- 0
poF[7,i,t,3] <- 0
poF[7,i,t,4] <- pDetT[t]*(1-Telem[i,t])+pDetT[t]*((Telem[i,t]*pFail[TelemYears[i,t],t]))
poF[7,i,t,5] <- 0
poF[7,i,t,6] <- 0        
poF[7,i,t,7] <- 0
poF[7,i,t,8] <- aAFC0O[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[7,i,t,9] <- 0
poF[7,i,t,10]<- 0
poF[7,i,t,11]<- aAFC0O[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[7,i,t,12]<- aAFC0O[3]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[7,i,t,13]<- 0
poF[7,i,t,14]<- 0
poF[7,i,t,15]<- 0
poF[7,i,t,16]<- 0
poF[7,i,t,17]<- (1-pDetT[t])*(1-Telem[i,t])+(1-pDetT[t])*((Telem[i,t]*pFail[TelemYears[i,t],t]))

poF[8,i,t,1] <- 0                  #AFC0/U
poF[8,i,t,2] <- 0
poF[8,i,t,3] <- 0
poF[8,i,t,4] <- 0
poF[8,i,t,5] <- 0
poF[8,i,t,6] <- 0
poF[8,i,t,7] <- 0
poF[8,i,t,8] <- 0
poF[8,i,t,9] <- 0
poF[8,i,t,10]<- 0
poF[8,i,t,11]<- aAFC0U[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[8,i,t,12]<- 0
poF[8,i,t,13]<- aAFC0U[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[8,i,t,14]<- 0
poF[8,i,t,15]<- 0
poF[8,i,t,16]<- 0
poF[8,i,t,17]<- (1-Telem[i,t])+(Telem[i,t]*pFail[TelemYears[i,t],t])

poF[9,i,t,1] <- 0                  #AFC1/O
poF[9,i,t,2] <- 0
poF[9,i,t,3] <- 0
poF[9,i,t,4] <- 0
poF[9,i,t,5] <- pDetT[t]*(1-Telem[i,t])+pDetT[t]*((Telem[i,t]*pFail[TelemYears[i,t],t]))
poF[9,i,t,6] <- 0
poF[9,i,t,7] <- 0
poF[9,i,t,8] <- 0
poF[9,i,t,9] <- aAFC1O[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[9,i,t,10]<- 0
poF[9,i,t,11]<- 0
poF[9,i,t,12]<- 0
poF[9,i,t,13]<- 0
poF[9,i,t,14]<- aAFC1O[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[9,i,t,15]<- aAFC1O[3]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[9,i,t,16]<- 0
poF[9,i,t,17]<- (1-pDetT[t])*(1-Telem[i,t])+(1-pDetT[t])*((Telem[i,t]*pFail[TelemYears[i,t],t]))

poF[10,i,t,1] <- 0                  #AFC1/U
poF[10,i,t,2] <- 0
poF[10,i,t,3] <- 0
poF[10,i,t,4] <- 0
poF[10,i,t,5] <- 0
poF[10,i,t,6] <- 0
poF[10,i,t,7] <- 0
poF[10,i,t,8] <- 0
poF[10,i,t,9] <- 0
poF[10,i,t,10]<- 0
poF[10,i,t,11]<- 0
poF[10,i,t,12]<- 0
poF[10,i,t,13]<- 0
poF[10,i,t,14]<- aAFC1U[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[10,i,t,15]<- 0
poF[10,i,t,16]<- aAFC1U[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[10,i,t,17]<- (1-Telem[i,t])+(Telem[i,t]*pFail[TelemYears[i,t],t])
        
poF[11,i,t,1] <- 0                  #AFC2/O
poF[11,i,t,2] <- 0
poF[11,i,t,3] <- 0
poF[11,i,t,4] <- 0
poF[11,i,t,5] <- 0
poF[11,i,t,6] <- pDetT[t]*(1-Telem[i,t])+pDetT[t]*((Telem[i,t]*pFail[TelemYears[i,t],t]))
poF[11,i,t,7] <- 0
poF[11,i,t,8] <- 0
poF[11,i,t,9] <- 0
poF[11,i,t,10]<- aAFC2O[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[11,i,t,11]<- 0
poF[11,i,t,12]<- 0
poF[11,i,t,13]<- 0
poF[11,i,t,14]<- aAFC2O[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[11,i,t,15]<- aAFC2O[3]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[11,i,t,16]<- 0
poF[11,i,t,17]<- (1-pDetT[t])*(1-Telem[i,t])+(1-pDetT[t])*((Telem[i,t]*pFail[TelemYears[i,t],t]))
        
poF[12,i,t,1] <- 0                  #AFC2/U
poF[12,i,t,2] <- 0
poF[12,i,t,3] <- 0
poF[12,i,t,4] <- 0
poF[12,i,t,5] <- 0
poF[12,i,t,6] <- 0
poF[12,i,t,7] <- 0
poF[12,i,t,8] <- 0
poF[12,i,t,9] <- 0
poF[12,i,t,10]<- 0
poF[12,i,t,11]<- 0
poF[12,i,t,12]<- 0
poF[12,i,t,13]<- 0
poF[12,i,t,14]<- aAFC2U[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[12,i,t,15]<- 0
poF[12,i,t,16]<- aAFC2U[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[12,i,t,17]<- (1-Telem[i,t])+(Telem[i,t]*pFail[TelemYears[i,t],t])

poF[13,i,t,1] <- 0                  #DEAD
poF[13,i,t,2] <- 0
poF[13,i,t,3] <- 0
poF[13,i,t,4] <- 0
poF[13,i,t,5] <- 0
poF[13,i,t,6] <- 0
poF[13,i,t,7] <- 0
poF[13,i,t,8] <- 0
poF[13,i,t,9] <- 0
poF[13,i,t,10]<- 0
poF[13,i,t,11]<- aDEAD[1]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[13,i,t,12]<- 0
poF[13,i,t,13]<- 0
poF[13,i,t,14]<- aDEAD[2]*Telem[i,t]*(1-pFail[TelemYears[i,t],t])
poF[13,i,t,15]<- 0
poF[13,i,t,16]<- 0
poF[13,i,t,17]<- (1-Telem[i,t])+(Telem[i,t]*pFail[TelemYears[i,t],t])

#End Female matrices
}}#i,t


## Component II. CMR MALE state and observation matrices
# Male State Matrix
for(i in 1:nMale){
 for(t in (firstMale[i]+1):nyrs){
  #State process (draw zM[t] given zM[t-1])
   zM[i,t] ~ dcat(psM[zM[i,t-1],])
  #Events (draw event[t] given zM[t])
   yM[i,t] ~ dcat(poM[zM[i,t], t-1,])
 }#t
}#i


psM[1,1] <- 0                         #M2Y/O 
psM[1,2] <- 0
psM[1,3] <- phiMSA*psiM[1]
psM[1,4] <- phiMSA*(1-psiM[1])
psM[1,5] <- 0
psM[1,6] <- 0
psM[1,7] <- 0
psM[1,8] <- 0
psM[1,9] <- (1-phiMSA)

psM[2,1] <- 0                         #M2Y/U
psM[2,2] <- 0
psM[2,3] <- phiMSA*(1-psiM[2])
psM[2,4] <- phiMSA*psiM[2]
psM[2,5] <- 0
psM[2,6] <- 0
psM[2,7] <- 0
psM[2,8] <- 0
psM[2,9] <- (1-phiMSA)

psM[3,1] <- 0                         #M3Y/O
psM[3,2] <- 0
psM[3,3] <- 0
psM[3,4] <- 0
psM[3,5] <- phiMSA*psiM[1]
psM[3,6] <- phiMSA*(1-psiM[1])
psM[3,7] <- 0
psM[3,8] <- 0
psM[3,9] <- (1-phiMSA)

psM[4,1] <- 0                         #M3Y/U
psM[4,2] <- 0
psM[4,3] <- 0
psM[4,4] <- 0
psM[4,5] <- phiMSA*(1-psiM[2])
psM[4,6] <- phiMSA*psiM[2]
psM[4,7] <- 0
psM[4,8] <- 0
psM[4,9] <- (1-phiMSA)

psM[5,1] <- 0                         #M4Y/O
psM[5,2] <- 0
psM[5,3] <- 0
psM[5,4] <- 0
psM[5,5] <- 0
psM[5,6] <- 0
psM[5,7] <- phiMSA*psiM[1]
psM[5,8] <- phiMSA*(1-psiM[1])
psM[5,9] <- (1-phiMSA)

psM[6,1] <- 0                         #M4Y/U
psM[6,2] <- 0
psM[6,3] <- 0
psM[6,4] <- 0
psM[6,5] <- 0
psM[6,6] <- 0
psM[6,7] <- phiMSA*(1-psiM[2])
psM[6,8] <- phiMSA*psiM[2]
psM[6,9] <- 1-phiMSA

psM[7,1] <- 0                         #AMY/O
psM[7,2] <- 0
psM[7,3] <- 0
psM[7,4] <- 0
psM[7,5] <- 0
psM[7,6] <- 0
psM[7,7] <- phiMAD*psiM[1]
psM[7,8] <- phiMAD*(1-psiM[1])
psM[7,9] <- 1-phiMAD

psM[8,1] <- 0                         #AMY/U
psM[8,2] <- 0
psM[8,3] <- 0
psM[8,4] <- 0
psM[8,5] <- 0
psM[8,6] <- 0
psM[8,7] <- phiMAD*(1-psiM[2])
psM[8,8] <- phiMAD*psiM[2]
psM[8,9] <- 1-phiMAD

psM[9,1] <- 0                         #DEAD
psM[9,2] <- 0
psM[9,3] <- 0
psM[9,4] <- 0
psM[9,5] <- 0
psM[9,6] <- 0
psM[9,7] <- 0
psM[9,8] <- 0
psM[9,9] <- 1


## Male Observation matrix
for(t in 1:(nyrs-1)){
poM[1,t,1] <- pDetT[t]                  #M2Y/O        
poM[1,t,2] <- 0
poM[1,t,3] <- 0
poM[1,t,4] <- 0
poM[1,t,5] <- (1-pDetT[t])

poM[2,t,1] <- 0                         #M2Y/U
poM[2,t,2] <- 0
poM[2,t,3] <- 0
poM[2,t,4] <- 0
poM[2,t,5] <- 1

poM[3,t,1] <- 0                         #M3Y/O
poM[3,t,2] <- pDetT[t]
poM[3,t,3] <- 0
poM[3,t,4] <- 0
poM[3,t,5] <- (1-pDetT[t])

poM[4,t,1] <- 0                         #M3Y/U
poM[4,t,2] <- 0
poM[4,t,3] <- 0
poM[4,t,4] <- 0
poM[4,t,5] <- 1

poM[5,t,1] <- 0                         #M4Y/O
poM[5,t,2] <- 0
poM[5,t,3] <- pDetT[t]
poM[5,t,4] <- 0
poM[5,t,5] <- (1-pDetT[t])

poM[6,t,1] <- 0                         #M4Y/U
poM[6,t,2] <- 0
poM[6,t,3] <- 0
poM[6,t,4] <- 0
poM[6,t,5] <- 1

poM[7,t,1] <- 0                         #AMY/O
poM[7,t,2] <- 0
poM[7,t,3] <- 0
poM[7,t,4] <- pDetT[t]
poM[7,t,5] <- (1-pDetT[t])

poM[8,t,1] <- 0                         #AMY/U
poM[8,t,2] <- 0
poM[8,t,3] <- 0
poM[8,t,4] <- 0
poM[8,t,5] <- 1

poM[9,t,1] <- 0                         #DEAD
poM[9,t,2] <- 0
poM[9,t,3] <- 0
poM[9,t,4] <- 0
poM[9,t,5] <- 1

}#t

#### Component III. Individual cub survival
 # The probability at least one C0 or C1 cub survived (phiL[1:2]), i.e. litter survival in cap-recap data,
 # is a function of individual C0 and C1 cub survial (phiCub[1:2])
 # and the proportion of 1,2,or 3 C0 and C1 cub litters (R0[1:3], R1[1:3])
  phiL0<-1-((1-phiCub0)^1*R0[1]+(1-phiCub0)^2*R0[2]+(1-phiCub0)^3*R0[3])
  phiL1<-1-((1-phiCub1)^1*R1[1]+(1-phiCub1)^2*R1[2]+(1-phiCub1)^3*R1[3])

#### Component IV. Litter size
 # The probability an AFC1 has 1,2, or 3 C1 cubs (R1[1:3,t])is a function of
 # The probability an AFC0 had 1,2, or 3 C0 cubs (R0[1:3]) and C0 cub survival (phiCub0[t])
 # NOTE: R1all[4,t] is the probability that entire litter died, which is already estimated above (1-phiL1[t]) 
  R1[1]<-(((phiCub0^1*R0[1])+(((1-phiCub0)*phiCub0*R0[2])*2)+((1-phiCub0)^2*phiCub0*R0[3])*3))/phiL0
  R1[2]<-(((phiCub0^2*R0[2])+(((1-phiCub0)*phiCub0^2*R0[3])*3)))/phiL0
  R1[3]<-((phiCub0^3*R0[3]))/phiL0

 # Observed number of C1 cubs (nObsC1yr) follows a categorical distribution
 for(t in 1:(nyrs-1)){
  nObsC1yr[1:nC1max,t]~dmulti(R1[],nObsAFC1yr[t]) 
 }#t

 # Dirichlet prior on R0[]
 for(k in 1:nC0max){
  R0p[k]~dgamma(1,1)
  R0[k]<-R0p[k]/sum(R0p[])
 }#k

 C0litterSize<-R0[1]+2*R0[2]+3*R0[3]
 C1litterSize<-R1[1]+2*R1[2]+3*R1[3]

#### Component V.  Additional data on weaning probability
# number of independent C2 (nWean) and total number of observed C2s (nC2obs) 
for(t in 1:nyrsW){
 nWean[t]~dbin(W, nC2obs[t]) 
}#t


### Component VI. Abundance and count data: Matrix-based abundance process model and state-specific count data
#state codes for abundance process
# 1 = F2Y/O
# 2 = F2Y/U
# 3 = F3Y/O
# 4 = F3Y/U
# 5 = AFNC/O
# 6 = AFNC/U
# 7 = AFC0/O
# 8 = AFC0/U
# 9 = AFC1/O
# 10 = AFC1/U
# 11 = AFC2/O
# 12 = AFC2/U
# 13 = M2Y/O
# 14 = M2Y/U
# 15 = M3Y/O
# 16 = M3Y/U
# 17 = M4Y/O
# 18 = M4Y/U
# 19 = AM/O
# 20 = AM/U


## Year 1 abundances for total stage (Ns_tot) and state-specific (Nstate) assuming SSD 
#Adult males are reference class
 Nstate[19,1] ~dpois(lambda[19])   #AM in         
 Nstate[20,1] ~dpois(lambda[20])   #AM out 
 ones[1]~dsum(step(Nstate[20,1]-Nstate[19,1]))   #force AM/out >= AM/In     
 Ns_tot[nClass]<-Nstate[19,1]+Nstate[20,1] 

#all other first year abundances are a function of AM and SSD 
for(s in c(1:(nClass-1))){ #using state AM as reference 
  Ns_tot[s]~dpois(Ns_tot[nClass]*ssd[s]/ssd[nClass])
  Nstate[(s*2-1),1] ~dpois(lambda[(s*2-1)])       #observable           
  Nstate[(s*2),1] <- Ns_tot[s]- Nstate[(s*2-1),1] #unobservable is derived   
} 

 #abundance priors for observable states and outside reference state (AM_out) in year 1
 for(s in c(1,3,5,7,9,11,13,15,17,19,20)){       
  lambda[s] ~ dgamma(0.001,0.001)             
 }

## Count data detection
for(t in 1){
 for(s in 1:nStages){
  nCounts[s,t] ~ dbin(pDet*ObsStage[s]*surveyYears[t], Nstate[s,t])
 }

 #Derive of dependent cubs
  #write multinomial as conditional binomials since JAGS does not allow latent N in dmulti
   #dependent C0 (sexes combined)
   AFC0O_1[t]~dbin(R0[1],Nstate[7,t])                        #adults with one age-0 cub
   AFC0O_2[t]~dbin(R0[2]/(1-R0[1]),(Nstate[7,t]-AFC0O_1[t])) #adults with two age-0 cubs
   AFC0O_3[t]<-Nstate[7,t]-AFC0O_1[t]-AFC0O_2[t]             #adults with three age-0 cubs
   Nc0O[t]<-round(AFC0O_1[t]+2*AFC0O_2[t]+3*AFC0O_3[t])      #age-0 cubs/in

   AFC0U_1[t]~dbin(R0[1],Nstate[8,t])                        #adults with one age-0 cub
   AFC0U_2[t]~dbin(R0[2]/(1-R0[1]),(Nstate[8,t]-AFC0U_1[t])) #adults with two age-0 cubs
   AFC0U_3[t]<-Nstate[8,t]-AFC0U_1[t]-AFC0U_2[t]             #adults with three age-0 cubs
   Nc0U[t]<-round(AFC0U_1[t]+2*AFC0U_2[t]+3*AFC0U_3[t] )     #age-0 cubs/in
  
   #dependent C1  (sexes combined)
   AFC1O_1[t]~dbin(R1[1],Nstate[9,t])                        #adults with one age-1 cub
   AFC1O_2[t]~dbin(R1[2]/(1-R1[1]),(Nstate[9,t]-AFC1O_1[t])) #adults with two age-1 cubs
   AFC1O_3[t]<-Nstate[9,t]-AFC1O_1[t]-AFC1O_2[t]             #adults with three age-1 cubs
   Nc1O[t]<-round(AFC1O_1[t]+2*AFC1O_2[t]+3*AFC1O_3[t])      #age-1 cubs/in

   AFC1U_1[t]~dbin(R1[1],Nstate[10,t])                        #adults with one age-1 cub
   AFC1U_2[t]~dbin(R1[2]/(1-R1[1]),(Nstate[10,t]-AFC1U_1[t])) #adults with two age-1 cubs
   AFC1U_3[t]<-Nstate[10,t]-AFC1U_1[t]-AFC1U_2[t]             #adults with three age-1 cubs
   Nc1U[t]<-round(AFC1U_1[t]+2*AFC1U_2[t]+3*AFC1U_3[t])       #age-1 cubs/in

 # total number of individuals exposed to sampling
 Nin[t]<-Nc0O[t]+Nc1O[t]+sum(Nstate[c(1,3,5,7,9,11,13,15,17,19),t]) #observable population
}

for(t in 2:nyrs){
 #count (detection) data 
 for(s in 1:nStages){
  nCounts[s,t] ~ dbin(pDet*ObsStage[s]*surveyYears[t], Nstate[s,t])
 }

 # total number of individuals exposed to sampling
 Nin[t]<-Nc0O[t]+Nc1O[t]+sum(Nstate[c(1,3,5,7,9,11,13,15,17,19),t]) #observable population


 #abundance process
  #dependent C0 (sexes combined)
  AFC0O_1[t]~dbin(R0[1],Nstate[7,t])                        #adults with one age-0 cub
  AFC0O_2[t]~dbin(R0[2]/(1-R0[1]),(Nstate[7,t]-AFC0O_1[t])) #adults with two age-0 cubs
  AFC0O_3[t]<-Nstate[7,t]-AFC0O_1[t]-AFC0O_2[t]             #adults with three age-0 cubs
  Nc0O[t]<-round(AFC0O_1[t]+2*AFC0O_2[t]+3*AFC0O_3[t])      #age-0 cubs/in

  AFC0U_1[t]~dbin(R0[1],Nstate[8,t])                        #adults with one age-0 cub
  AFC0U_2[t]~dbin(R0[2]/(1-R0[1]),(Nstate[8,t]-AFC0U_1[t])) #adults with two age-0 cubs
  AFC0U_3[t]<-Nstate[8,t]-AFC0U_1[t]-AFC0U_2[t]             #adults with three age-0 cubs
  Nc0U[t]<-round(AFC0U_1[t]+2*AFC0U_2[t]+3*AFC0U_3[t] )     #age-0 cubs/in
  
  #dependent C1  (sexes combined)
  AFC1O_1[t]~dbin(R1[1],Nstate[9,t])                        #adults with one age-1 cub
  AFC1O_2[t]~dbin(R1[2]/(1-R1[1]),(Nstate[9,t]-AFC1O_1[t])) #adults with two age-1 cubs
  AFC1O_3[t]<-Nstate[9,t]-AFC1O_1[t]-AFC1O_2[t]             #adults with three age-1 cubs
  Nc1O[t]<-round(AFC1O_1[t]+2*AFC1O_2[t]+3*AFC1O_3[t])      #age-1 cubs/in

  AFC1U_1[t]~dbin(R1[1],Nstate[10,t])                        #adults with one age-1 cub
  AFC1U_2[t]~dbin(R1[2]/(1-R1[1]),(Nstate[10,t]-AFC1U_1[t])) #adults with two age-1 cubs
  AFC1U_3[t]<-Nstate[10,t]-AFC1U_1[t]-AFC1U_2[t]             #adults with three age-1 cubs
  Nc1U[t]<-round(AFC1U_1[t]+2*AFC1U_2[t]+3*AFC1U_3[t])       #age-1 cubs/in

  #C1O transition process 
   Nc1O_2YO[t-1]~dbin(phiFAD*phiCub1*psiF[1],     Nc1O[t-1])   #age 1/in to 2YO/in (sexes combined)
   Nc1O_2YU[t-1]~dbin(phiFAD*phiCub1*(1-psiF[1]), Nc1O[t-1])   #age 1/in to 2YO/out (sexes combined)
   Nc1O_Dead[t-1]<-Nc1O[t-1]-(Nc1O_2YO[t-1]+Nc1O_2YU[t-1])     #age 1/in to Dead
   onesC1[1,t-1]~dsum(step(Nc1O_Dead[t-1]))                    #number dead cannot be negative

  #C1U transition process 
   Nc1U_2YO[t-1]~dbin(phiFAD*phiCub1*(1-psiF[2]),Nc1U[t-1])    #age 1/out to 2YO/in (sexes combined)
   Nc1U_2YU[t-1]~dbin(phiFAD*phiCub1*psiF[2],    Nc1U[t-1])    #age 1/out to 2YO/out (sexes combined)
   Nc1U_Dead[t-1]<-Nc1U[t-1]-(Nc1U_2YO[t-1]+Nc1U_2YU[t-1])     #age 1/in to Dead
   onesC1[2,t-1]~dsum(step(Nc1U_Dead[t-1]))                    #number dead cannot be negative

### FEMALE MATRICES
 #F2Y
 Nstate[1,t]~dbin(.5,Nc1O_2YO[t-1] + Nc1U_2YO[t-1])           #F2Y/O (.5 for equal sex ratio)
 Nstate[2,t]~dbin(.5,Nc1O_2YU[t-1] + Nc1U_2YU[t-1])           #F2Y/U (.5 for equal sex ratio)

  #F2Y transitions
   F2YO_F3YO[t-1]~dbin(phiFSA*psiF[1], Nstate[1,t-1])           #F2Y/in to F2Y/in
   F2YO_F3YU[t-1]~dbin(phiFSA*(1-psiF[1]), Nstate[1,t-1])       #F2Y/in to F2Y/out 
   Ndead[1,t-1]<-Nstate[1,t-1]-(F2YO_F3YO[t-1]+F2YO_F3YU[t-1])  #F2Y/in to dead 
   onesN[1,t-1]~dsum(step(Ndead[1,t-1]))                        #number dead cannot be negative

   F2YU_F3YO[t-1]~dbin(phiFSA*(1-psiF[2]), Nstate[2,t-1])          #F2Y/out to F2Y/in
   F2YU_F3YU[t-1]~dbin(phiFSA*psiF[2],     Nstate[2,t-1])          #F2Y/out to F2Y/out
   Ndead[2,t-1]<-Nstate[2,t-1]-(F2YU_F3YO[t-1]+F2YU_F3YU[t-1])     #F2Y/out to dead
   onesN[2,t-1]~dsum(step(Ndead[2,t-1]))                           #number dead cannot be negative

 #F3Y
 Nstate[3,t] <- F2YO_F3YO[t-1] + F2YU_F3YO[t-1]        #F2Y/in that survived and stayed in + F2Y/out that survived and moved in  
 Nstate[4,t] <- F2YO_F3YU[t-1] + F2YU_F3YU[t-1]        #F2Y/in that survived and moved out + F2Y/out that survived and stayed out  

  #F3Y transitions
   F3YO_AFNCO[t-1]~dbin(phiFSA*psiF[1],     Nstate[3,t-1])       #F3Y/in to AFNC/in
   F3YO_AFNCU[t-1]~dbin(phiFSA*(1-psiF[1]), Nstate[3,t-1])       #F3Y/in to AFNC/out
   Ndead[3,t-1]<-Nstate[3,t-1]-(F3YO_AFNCO[t-1]+F3YO_AFNCU[t-1]) #number dead 
   onesN[3,t-1]~dsum(step(Ndead[3,t-1]))                         #number dead cannot be negative

   F3YU_AFNCO[t-1]~dbin(phiFSA*(1-psiF[2]), Nstate[4,t-1])       #F3Y/out to AFNC/in
   F3YU_AFNCU[t-1]~dbin(phiFSA*psiF[2],     Nstate[4,t-1])       #F3Y/out to AFNC/out
   Ndead[4,t-1]<-Nstate[4,t-1]-(F3YU_AFNCO[t-1]+F3YU_AFNCU[t-1]) #number dead
   onesN[4,t-1]~dsum(step(Ndead[4,t-1]))                         #number dead cannot be negative
  
 #AFNC (this gets ridiculous due to all the transition possibilities)
 Nstate[5,t] <- F3YO_AFNCO[t-1] + F3YU_AFNCO[t-1] + AFNCO_AFNCO[t-1] + AFNCU_AFNCO[t-1] + AFC0O_AFNCO[t-1] + AFC0U_AFNCO[t-1] + AFC1O_AFNCO[t-1] + AFC1U_AFNCO[t-1] + AFC2O_AFNCO[t-1] + AFC2U_AFNCO[t-1]
 Nstate[6,t] <- F3YO_AFNCU[t-1] + F3YU_AFNCU[t-1] + AFNCO_AFNCU[t-1] + AFNCU_AFNCU[t-1] + AFC0O_AFNCU[t-1] + AFC0U_AFNCU[t-1] + AFC1O_AFNCU[t-1] + AFC1U_AFNCU[t-1] + AFC2O_AFNCU[t-1] + AFC2U_AFNCU[t-1]

  #AFNCO transitions
   AFNCO_AFNCO[t-1]~dbin(phiFAD*(1-B1)*psiF[1],     Nstate[5,t-1])  #AFNC/in to AFNC/in
   AFNCO_AFNCU[t-1]~dbin(phiFAD*(1-B1)*(1-psiF[1]), Nstate[5,t-1])  #AFNC/in to AFNC/out
   AFNCO_AFC0O[t-1]~dbin(phiFAD*B1*psiF[3],         Nstate[5,t-1])  #AFNC/in to AFC0/in
   AFNCO_AFC0U[t-1]~dbin(phiFAD*B1*(1-psiF[3]),     Nstate[5,t-1])  #AFNC/in to AFC0/in
   Ndead[5,t-1]<-Nstate[5,t-1]-(AFNCO_AFNCO[t-1]+AFNCO_AFNCU[t-1]+AFNCO_AFC0O[t-1]+AFNCO_AFC0U[t-1])  #number dead 
   onesN[5,t-1]~dsum(step(Ndead[5,t-1]))                            #number dead cannot be negative

  #AFNCU transitions
   AFNCU_AFNCO[t-1]~dbin(phiFAD*(1-B1)*(1-psiF[2]), Nstate[6,t-1])    #AFNC/out to AFNC/in
   AFNCU_AFNCU[t-1]~dbin(phiFAD*(1-B1)*psiF[2],     Nstate[6,t-1])    #AFNC/out to AFNC/out
   AFNCU_AFC0O[t-1]~dbin(phiFAD*B1*(1-psiF[4]),     Nstate[6,t-1])    #AFNC/out to AFC0/in
   AFNCU_AFC0U[t-1]~dbin(phiFAD*B1*psiF[4],         Nstate[6,t-1])    #AFNC/out to AFC0/out
   Ndead[6,t-1]<-Nstate[6,t-1]-( AFNCU_AFNCO[t-1]+AFNCU_AFNCU[t-1]+AFNCU_AFC0O[t-1]+AFNCU_AFC0U[t-1])  #number dead 
   onesN[6,t-1]~dsum(step(Ndead[6,t-1]))                              #number dead cannot be negative

 #AFC0
 Nstate[7,t] <- AFNCO_AFC0O[t-1] + AFNCU_AFC0O[t-1] + AFC0O_AFC0O[t-1] + AFC0U_AFC0O[t-1] + AFC1O_AFC0O[t-1] + AFC1U_AFC0O[t-1] + AFC2O_AFC0O[t-1] + AFC2U_AFC0O[t-1] 
 Nstate[8,t] <- AFNCO_AFC0U[t-1] + AFNCU_AFC0U[t-1] + AFC0O_AFC0U[t-1] + AFC0U_AFC0U[t-1] + AFC1O_AFC0U[t-1] + AFC1U_AFC0U[t-1] + AFC2O_AFC0U[t-1] + AFC2U_AFC0U[t-1]

  #AFC0O transitions
   AFC0O_AFC1O[t-1]~dbin(phiFAD*phiL0*psiF[1],                     Nstate[7,t-1])    #AFC0/in to AFC1/in
   AFC0O_AFC1U[t-1]~dbin(phiFAD*phiL0*(1-psiF[1]),                 Nstate[7,t-1])    #AFC0/in to AFC1/out
   AFC0O_AFC0O[t-1]~dbin(phiFAD*(1-phiL0)*B2*psiF[3],              Nstate[7,t-1])    #AFC0/in to AFC0/in
   AFC0O_AFC0U[t-1]~dbin(phiFAD*(1-phiL0)*B2*(1-psiF[3]),          Nstate[7,t-1])    #AFC0/in to AFC0/out
   AFC0O_AFNCO[t-1]~dbin(phiFAD*(1-phiL0)*(1-B2)*psiF[1],          Nstate[7,t-1])    #AFC0/in to AFNC/in
   AFC0O_AFNCU[t-1]~dbin(phiFAD*(1-phiL0)*(1-B2)*(1-psiF[1]),      Nstate[7,t-1])    #AFC0/in to AFNC/out
   Ndead[7,t-1]<-Nstate[7,t-1]-(AFC0O_AFC1O[t-1]+AFC0O_AFC1U[t-1]+AFC0O_AFC0O[t-1]+AFC0O_AFC0U[t-1]+AFC0O_AFNCO[t-1]+AFC0O_AFNCU[t-1])  #number dead 
   onesN[7,t-1]~dsum(step(Ndead[7,t-1]))                                            #number dead cannot be negative

  #AFC0U transitions
   AFC0U_AFC1O[t-1]~dbin(phiFAD*phiL0*(1-psiF[2]),                 Nstate[8,t-1])    #AFC0/out to AFC1/in
   AFC0U_AFC1U[t-1]~dbin(phiFAD*phiL0*psiF[2],                     Nstate[8,t-1])    #AFC0/out to AFC1/out
   AFC0U_AFC0O[t-1]~dbin(phiFAD*(1-phiL0)*B2*(1-psiF[4]),          Nstate[8,t-1])    #AFC0/out to AFC0/in
   AFC0U_AFC0U[t-1]~dbin(phiFAD*(1-phiL0)*B2*psiF[4],              Nstate[8,t-1])    #AFC0/out to AFC0/out
   AFC0U_AFNCO[t-1]~dbin(phiFAD*(1-phiL0)*(1-B2)*(1-psiF[2]),      Nstate[8,t-1])    #AFC0/out to AFNC/in
   AFC0U_AFNCU[t-1]~dbin(phiFAD*(1-phiL0)*(1-B2)*psiF[2],          Nstate[8,t-1])    #AFC0/out to AFNC/out
   Ndead[8,t-1]<-Nstate[8,t-1]-(AFC0U_AFC1O[t-1]+AFC0U_AFC1U[t-1]+AFC0U_AFC0O[t-1]+AFC0U_AFC0U[t-1]+AFC0U_AFNCO[t-1]+AFC0U_AFNCU[t-1])  #number dead 
   onesN[8,t-1]~dsum(step(Ndead[8,t-1]))                                            #number dead cannot be negative

 #AFC1
 Nstate[9,t] <- AFC0O_AFC1O[t-1] + AFC0U_AFC1O[t-1] 
 Nstate[10,t] <- AFC0O_AFC1U[t-1] + AFC0U_AFC1U[t-1]

  #AFC1O transitions
   AFC1O_AFC2O[t-1]~dbin(phiFAD*phiL1*(1-W)*psiF[1],                         Nstate[9,t-1])       #AFC1/in to AFC2/in
   AFC1O_AFC2U[t-1]~dbin(phiFAD*phiL1*(1-W)*(1-psiF[1]),                     Nstate[9,t-1])       #AFC1/in to AFC2/out
   AFC1O_AFC0O[t-1]~dbin(phiFAD*psiF[3]*B2*(1-phiL1),                        Nstate[9,t-1])       #AFC1/in to AFC0/in
   AFC1O_AFC0U[t-1]~dbin(phiFAD*(1-psiF[3])*B2*(1-phiL1),                    Nstate[9,t-1])       #AFC1/in to AFC0/out
   AFC1O_AFNCO[t-1]~dbin(phiFAD*psiF[1]*((phiL1*W)+((1-B2)*(1-phiL1))),      Nstate[9,t-1])       #AFC1/in to AFNC/in
   AFC1O_AFNCU[t-1]~dbin(phiFAD*(1-psiF[1])*((phiL1*W)+(1-B2)*(1-phiL1)),    Nstate[9,t-1])       #AFC1/in to AFNC/out
   Ndead[9,t-1]<-Nstate[9,t-1]-(AFC1O_AFC2O[t-1]+AFC1O_AFC2U[t-1]+AFC1O_AFC0O[t-1]+AFC1O_AFC0U[t-1]+AFC1O_AFNCO[t-1]+AFC1O_AFNCU[t-1])  #number dead 
   onesN[9,t-1]~dsum(step(Ndead[9,t-1]))                                                         #number dead cannot be negative

  #AFC1U transitions
   AFC1U_AFC2O[t-1]~dbin(phiFAD*phiL1*(1-W)*(1-psiF[2]),                     Nstate[10,t-1])       #AFC1/out to AFC2/in
   AFC1U_AFC2U[t-1]~dbin(phiFAD*phiL1*(1-W)*psiF[2],                         Nstate[10,t-1])       #AFC1/out to AFC2/out
   AFC1U_AFC0O[t-1]~dbin(phiFAD*(1-psiF[4])*B2*(1-phiL1),                    Nstate[10,t-1])       #AFC1/out to AFC0/in
   AFC1U_AFC0U[t-1]~dbin(phiFAD*psiF[4]*B2*(1-phiL1),                        Nstate[10,t-1])       #AFC1/out to AFC0/out
   AFC1U_AFNCO[t-1]~dbin(phiFAD*(1-psiF[2])*((phiL1*W)+(1-B2)*(1-phiL1)),    Nstate[10,t-1])       #AFC1/out to AFNC/in
   AFC1U_AFNCU[t-1]~dbin(phiFAD*psiF[2]*((phiL1*W)+(1-B2)*(1-phiL1)),        Nstate[10,t-1])       #AFC1/out to AFNC/out
   Ndead[10,t-1]<-Nstate[10,t-1]-(AFC1U_AFC2O[t-1]+AFC1U_AFC2U[t-1]+AFC1U_AFC0O[t-1]+AFC1U_AFC0U[t-1]+AFC1U_AFNCO[t-1]+AFC1U_AFNCU[t-1])  #number dead 
   onesN[10,t-1]~dsum(step(Ndead[10,t-1]))                                                         #number dead cannot be negative

 #AFC2
 Nstate[11,t] <- AFC1O_AFC2O[t-1] + AFC1U_AFC2O[t-1]
 Nstate[12,t] <- AFC1O_AFC2U[t-1] + AFC1U_AFC2U[t-1]

  #AFC2O transitions
   AFC2O_AFC0O[t-1]~dbin(phiFAD*B1*psiF[3],          Nstate[11,t-1])                #AFC2/in to AFC0/in
   AFC2O_AFC0U[t-1]~dbin(phiFAD*B1*(1-psiF[3]),      Nstate[11,t-1])                #AFC2/in to AFC0/out
   AFC2O_AFNCO[t-1]~dbin(phiFAD*(1-B1)*psiF[1],      Nstate[11,t-1])                #AFC2/in to AFNC/in
   AFC2O_AFNCU[t-1]~dbin(phiFAD*(1-B1)*(1-psiF[1]),  Nstate[11,t-1])                #AFC2/in to AFNC/out
   Ndead[11,t-1]<-Nstate[11,t-1]-(AFC2O_AFC0O[t-1]+AFC2O_AFC0U[t-1]+AFC2O_AFNCO[t-1]+AFC2O_AFNCU[t-1])  #number dead 
   onesN[11,t-1]~dsum(step(Ndead[11,t-1]))  

  #AFC2U transitions
   AFC2U_AFC0O[t-1]~dbin(phiFAD*B1*(1-psiF[4]),      Nstate[12,t-1])                 #AFC2/out to AFC0/in
   AFC2U_AFC0U[t-1]~dbin(phiFAD*B1*psiF[4],          Nstate[12,t-1])                 #AFC2/out to AFC0/out
   AFC2U_AFNCO[t-1]~dbin(phiFAD*(1-B1)*(1-psiF[2]),  Nstate[12,t-1])                 #AFC2/out to AFNC/in
   AFC2U_AFNCU[t-1]~dbin(phiFAD*(1-B1)*psiF[2],      Nstate[12,t-1])                 #AFC2/out to AFNC/out
   Ndead[12,t-1]<-Nstate[12,t-1]-(AFC2U_AFC0O[t-1]+AFC2U_AFC0U[t-1]+AFC2U_AFNCO[t-1]+AFC2U_AFNCU[t-1])  #number dead 
   onesN[12,t-1]~dsum(step(Ndead[12,t-1]))  

### End Females


### MALE MATRICES
 #M2Y
 Nstate[13,t]<- (Nc1O_2YO[t-1] + Nc1U_2YO[t-1])-Nstate[1,t]       #M2Y/O are 2-yr olds/in that were not females (Nstate[1,t]) plus new immigrants
 Nstate[14,t]<- (Nc1O_2YU[t-1] + Nc1U_2YU[t-1])-Nstate[2,t]       #M2Y/U are 2-yr olds/out that were not females (Nstate[2,t])

 #M2Y/in transitions
   M2YO_M3YO[t-1]~dbin(phiMSA*psiM[1],     Nstate[13,t-1])        #M2Y/in to M3Y/in
   M2YO_M3YU[t-1]~dbin(phiMSA*(1-psiM[1]), Nstate[13,t-1])        #M2Y/in to M3Y/out
   Ndead[13,t-1]<-Nstate[13,t-1]-(M2YO_M3YO[t-1]+M2YO_M3YU[t-1])  #number dead 
   onesN[13,t-1]~dsum(step(Ndead[13,t-1]))  

  #M2Y/out transitions
   M2YU_M3YO[t-1]~dbin(phiMSA*(1-psiM[2]), Nstate[14,t-1])        #M2Y/out to M3Y/in
   M2YU_M3YU[t-1]~dbin(phiMSA*psiM[2],     Nstate[14,t-1])        #M2Y/out to M3Y/out
   Ndead[14,t-1]<-Nstate[14,t-1]-(M2YU_M3YO[t-1]+M2YU_M3YU[t-1])  #number dead 
   onesN[14,t-1]~dsum(step(Ndead[14,t-1]))  


 #M3Y
 Nstate[15,t]<- M2YO_M3YO[t-1] + M2YU_M3YO[t-1]                #M3Y/O 
 Nstate[16,t]<- M2YO_M3YU[t-1] + M2YU_M3YU[t-1]                #M3Y/U 

  #M3Y/in transitions
   M3YO_M4YO[t-1]~dbin(phiMSA*psiM[1],     Nstate[15,t-1])        #M3Y/in to M4Y/in
   M3YO_M4YU[t-1]~dbin(phiMSA*(1-psiM[1]), Nstate[15,t-1])        #M3Y/in to M4Y/out
   Ndead[15,t-1]<-Nstate[15,t-1]-(M3YO_M4YO[t-1]+M3YO_M4YU[t-1])  #number dead 
   onesN[15,t-1]~dsum(step(Ndead[15,t-1]))  

  #M3Y/out transitions
   M3YU_M4YO[t-1]~dbin(phiMSA*(1-psiM[2]), Nstate[16,t-1])     #M3Y/out to M4Y/in
   M3YU_M4YU[t-1]~dbin(phiMSA*psiM[2],     Nstate[16,t-1])     #M3Y/out to M4Y/out
   Ndead[16,t-1]<-Nstate[16,t-1]-(M3YU_M4YO[t-1]+M3YU_M4YU[t-1])  #number dead 
   onesN[16,t-1]~dsum(step(Ndead[16,t-1]))  


 #M4Y
 Nstate[17,t]<- M3YO_M4YO[t-1] + M3YU_M4YO[t-1]                #M4Y/O 
 Nstate[18,t]<- M3YO_M4YU[t-1] + M3YU_M4YU[t-1]                #M4Y/U 

  #M4Y/in transitions
   M4YO_AMO[t-1]~dbin(phiMSA*psiM[1],     Nstate[17,t-1])       #M4Y/in to AM/in
   M4YO_AMU[t-1]~dbin(phiMSA*(1-psiM[1]), Nstate[17,t-1])       #M4Y/in to AM/out
   Ndead[17,t-1]<-Nstate[17,t-1]-(M4YO_AMO[t-1]+M4YO_AMU[t-1])  #number dead 
   onesN[17,t-1]~dsum(step(Ndead[17,t-1]))  

  #M4Y/out transitions
   M4YU_AMO[t-1]~dbin(phiMSA*(1-psiM[2]), Nstate[18,t-1])       #M4Y/out to AM/in
   M4YU_AMU[t-1]~dbin(phiMSA*psiM[2],     Nstate[18,t-1])       #M4Y/out to AM/out
   Ndead[18,t-1]<-Nstate[18,t-1]-(M4YU_AMO[t-1]+M4YU_AMU[t-1])  #number dead 
   onesN[18,t-1]~dsum(step(Ndead[18,t-1]))  


 #Adult Male
 Nstate[19,t]<- M4YO_AMO[t-1] + M4YU_AMO[t-1] + AMO_AMO[t-1] + AMU_AMO[t-1]  #AdultMale/O 
 Nstate[20,t]<- M4YO_AMU[t-1] + M4YU_AMU[t-1] + AMO_AMU[t-1] + AMU_AMU[t-1]  #AdM/U 
 ones[t]~dsum(step(Nstate[20,t]-Nstate[19,t]))                               #force AM/out >= AM/In     

  #AM/in transitions
   AMO_AMO[t-1]~dbin(phiMAD*psiM[1],     Nstate[19,t-1])              #AM/in to AM/in
   AMO_AMU[t-1]~dbin(phiMAD*(1-psiM[1]), Nstate[19,t-1])              #AM/in to AM/out
   Ndead[19,t-1]<-Nstate[19,t-1]-(AMO_AMO[t-1]+AMO_AMU[t-1])          #number dead 
   onesN[19,t-1]~dsum(step(Ndead[19,t-1]))  

  #AdM/out transitions
   AMU_AMO[t-1]~dbin(phiMAD*(1-psiM[2]), Nstate[20,t-1])               #AM/out to AM/in
   AMU_AMU[t-1]~dbin(phiMAD*psiM[2],     Nstate[20,t-1])               #AM/out to AM/out
   Ndead[20,t-1]<-Nstate[20,t-1]-(AMU_AMO[t-1]+AMU_AMU[t-1])           #number dead 
   onesN[20,t-1]~dsum(step(Ndead[20,t-1]))   

}#nyrs

### End abundance process

#### Component VII.  Restrictions and Priors
# pFail[j = yr since deployment, t = release yr]
# pFail may vary by year 
for(t in 1:(nyrs-1)){
 pFail[1,t]<-0           # Tags do not fail prior to deployment year (j=1)
 pFail[2,t]~dunif(0,1)   # Tags may fail during first year of deployment (j=2)
 for(j in 3:MaxTelYr){
  pFail[j,t]<-1          # Tags are not active after first year of deployment(j>2)  
 }
}

#Survival priors
phiCub0~dunif(0,1) 
phiCub1~dunif(0,1)
phiFSA~dbeta(33.96, 4.20)   
phiFAD~dbeta(150.43, 11.32) 
phiMSA~dbeta(11.28, 2.48)   
phiMAD~dbeta(33.96, 4.20)   

#Breeding prob. priors and covariates
B1~dunif(0,1)
B2~dbeta(2.1,11.4)

#other priors
 W~dunif(0,1)         # weaning probability
 psiM[1]<-psiF[1]     # male transitions equal to female (not C0) transitions
 psiM[2]<-psiF[2]     # male transitions equal to female (not C0) transitions
 for (s in 1:4){
  psiF[s]~dunif(0,1)  
 }#s                  

#observation priors
  pDetT[1]<-pDet    #CMR det prob 2009
  pDetT[2]<-pDet    #CMR det prob 2010
  pDetT[3]<-pDet    #CMR det prob 2011
  pDetT[4]<-0       #CMR det prob 2012(no surveys)
  pDetT[5]<-pDet    #CMR det prob 2013
  pDetT[6]<-0       #CMR det prob 2014(no surveys)
  pDetT[7]<-pDet    #CMR det prob 2015
  pDetT[8]<-pDet    #CMR det prob 2016
  pDet~dunif(0,1)   #det prob constant across years

 #vague dirichlet priors for multi-event observation probabilities
 #indexing is messy since each telemetry event maps to different states
 for(i in 1:5){
   a5[i]~dgamma(1,1)
   aAFNCO[i]<-a5[i]/sum(a5[])
 }#i

 for(i in 1:4){
   a6[i]~dgamma(1,1)
   aAFNCU[i]<-a6[i]/sum(a6[])
 }#i

 for(i in 1:3){
   a7[i]~dgamma(1,1)
   aAFC0O[i]<-a7[i]/sum(a7[])
   a9[i]~dgamma(1,1)
   aAFC1O[i]<-a9[i]/sum(a9[])
   a11[i]~dgamma(1,1)
   aAFC2O[i]<-a11[i]/sum(a11[])
 }#i

 for(i in 1:2){
   a8[i]~dgamma(1,1)
   aAFC0U[i]<-a8[i]/sum(a8[])
   a10[i]~dgamma(1,1)
   aAFC1U[i]<-a10[i]/sum(a10[])
   a12[i]~dgamma(1,1)
   aAFC2U[i]<-a12[i]/sum(a12[])
   a13[i]~dgamma(1,1)
   aDEAD[i]<-a13[i]/sum(a13[])
 }#i



## Fit for count data
for(t in 1:nyrs){
 for(s in 1:nStages){
  eval[s,t] <- pDet*ObsStage[s]*surveyYears[t] * Nstate[s,t]     # Expected values
  y.new[s,t] ~ dbin(pDet*ObsStage[s]*surveyYears[t], Nstate[s,t]) # Generate replicate data

  # Chi-squared discrepancy with small constant to prevent division by zero
  CS_E[s,t] <-     pow((nCounts[s,t] - eval[s,t]),2) / (eval[s,t] + 0.5)
  CS_E.new[s,t] <- pow((y.new[s,t] - eval[s,t]),2) / (eval[s,t] + 0.5)

 }
}

CS_N_fit <- sum(CS_E[,])
CS_N_fitNew <- sum(CS_E.new[,])


## Fit for litter size data
for(t in 1:(nyrs-1)){
  LS_new[1:3,t] ~ dmulti(R1[1:3]*surveyYears[t+1],nObsAFC1yr[t]) # Generate replicate data
  for(s in 1:3){
   LS_eval[s,t] <- R1[s]*surveyYears[t+1] * nObsAFC1yr[t]   # Expected values

   #Chi-squared
   CS_LS_E[s,t] <- pow((nObsC1yr[s,t] - LS_eval[s,t]),2) / (LS_eval[s,t] + 0.5)
   CS_LS_E.new[s,t] <- pow((LS_new[s,t] - LS_eval[s,t]),2) / (LS_eval[s,t] + 0.5)

   }#s
}#t

CS_LS_fit <- sum(CS_LS_E[,])
CS_LS_fitNew <- sum(CS_LS_E.new[,])


## Fit for weaning data
for(t in 1:nyrsW){
  NW_eval[t] <-  W*nC2obs[t]      # Expected values
  NW_new[t] ~ dbin(W, nC2obs[t])  # Generate replicated data

  #Chi-squared
  CS_NW_E[t] <- pow((nWean[t] - NW_eval[t]),2) / (NW_eval[t] + 0.5)
  CS_NW_E.new[t] <- pow((NW_new[t] - NW_eval[t]),2) / (NW_eval[t] + 0.5)
}

CS_NW_fit <-    sum(CS_NW_E[])
CS_NW_fitNew <- sum(CS_NW_E.new[])


}#model
",file = "Regehr2018_PolarBear_IPM_JAGS_model_GOF.txt")



##
##
##
## END JAGS MODEL
##
##
##



