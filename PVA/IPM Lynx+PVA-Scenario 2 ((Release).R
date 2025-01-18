#==============================================================================#
#                                                                              #
#                    IBERIAN LYNX INTEGRATED POPULATION MODEL                  #
#                            ~~~ Lynx pardinus ~~~                             #
#                                PVA-Releases                                  #
#       José Jiménez, Matías Taborda, Pablo Ferreras, Maria Jesús Palacios,    #
#        Fernando David Nájera, Jorge Peña, Marc Kéry and Michael Schaub       # 
#                            09/08/2023 9:41:42                                #
#                                                                              #
#==============================================================================#

# packages
library(IPMbook)
library(MCMCvis)
library(nimble)
library(mcmcOutput)

setwd('...')

bdataH <- read.csv('females.csv', header=TRUE, sep=";")
bdataH<-bdataH[bdataH$Tras=='0',]

# Released individuals (1: wild-born; 2: released)
suelta1 <- as.numeric(bdataH$suelta)
CH1<-data.matrix(bdataH[,9:19])

# Replace the 12 with a 0, because 12 means 'not seen'
CH1[CH1==12] <- 0

# Compute the multistate m-array from the individual capture histories, including one unobservable state
marrF <- marray(CH1, unobs=1, groups=suelta1)


bdataM <- read.csv('males.csv', header=TRUE, sep=";")
bdataM<-bdataM[bdataM$Tras=='0',]

# Released individuals (1: wild-born; 2: released)
suelta2 <- as.numeric(bdataM$suelta)
CH2<-data.matrix(bdataM[,9:19])

# Replace the 10 with a 0, because 10 means 'not seen'
CH2[CH2==10] <- 0

marrM <- marray(CH2, unobs=1, groups=suelta2)

# We first need to define two nimble functions that help with the calculation of the multistate m-array cell probabilities

###############################
#
# Function to perform matrix multiplication of multiple matrices
#
# Input:
# - mat: Array with dimension = {nstates, nstates, q}
# 
# - Note: q is the number of matrices that are multiplied in order: 1 * 2 * ... * q
#
# Output:
# - out: Matrix with dimension = {nstates, nstates}
#
# Created: 16.3.2023
# Author: Michael Schaub
#
##############################


matrixMult <- nimbleFunction(
  run = function(mat=double(3)) {
  di <- dim(mat)
  out <- matrix(NA, ncol=di[1], nrow=di[1])
  out <- mat[,,1] %*% mat[,,2]
  if (nimDim(mat)[3] > 2) {
    for (t in 3:di[3]){
      out <- out %*% mat[,,t]
    }
  }
  return(out)
  returnType(double(2))
})


###############################
#
# Function to calculate the probabilities of a multistate m-array based on transition and recapture probabilities
#
# Input:
# - psi: transition probabilities. Array with dimension = {nstates, nyears-1, nstates}
# - dp: recapture probabilities. Array with dimension = {nstates, nyears-1, nstates}
# - dq: probabilities of not being recaptured. Array with dimension = {nstates, nyears-1, nstates}
#
# Output:
# - pi: probabilities of the multistate m-array. Matrix with dimension = {nstates * (nyears-1), nyears*nstates-(nstates-1)}
#
# Created: 16.3.2023
# Author: Michael Schaub
#
##############################

marrProb <- nimbleFunction(
  run = function(psi=double(3), dp=double(3), dq=double(3)) {
  di <- dim(psi)
  ns <- di[1]                     # number of states
  nyears <- di[3] + 1             # number of years (occasions)
  pi <- array(0, dim=c((nyears-1)*ns, nyears*ns-(ns-1)))    # marray probabilities
  index <- matrix(NA, nrow=nyears-1, ncol=2)                # indices for m-array position for each occasion
  for (t in 1:(nyears-1)){
    index[t,1] <- (t-1)*ns + 1
    index[t,2] <- (t-1)*ns + ns
  }
  
  # Calculate cell probabilities
  # Diagonal
  for (t in 1:(nyears-1)){
    pi[index[t,1]:index[t,2], index[t,1]:index[t,2]] <- psi[1:ns,1:ns,t] %*% dp[1:ns,1:ns,t]	
  } #t
  # Above main diagonal
  for (t in 1:(nyears-2)){
    arr <- array(NA, dim=c(ns,ns,(nyears-1)*2))       # contains the transition and recpature matrices to be multiplied
    for (k in 1:((nyears-1)*2)){
      arr[1:ns,1:ns,k] <- diag(ns)
    } #k
    for (j in (t+1):(nyears-1)){
      for (m in t:(j-1)){
        arr[1:ns,1:ns,(2*m-1)] <- psi[1:ns,1:ns,m]
        arr[1:ns,1:ns,2*m] <- dq[1:ns,1:ns,m]
      } #m
    arr[1:ns,1:ns,(2*j-1)] <- psi[1:ns,1:ns,j]
    arr[1:ns,1:ns,2*j] <- dp[1:ns,1:ns,j]
    pi[index[t,1]:index[t,2], index[j,1]:index[j,2]] <- matrixMult(arr[,,])
    } #j
  } #t
  # Probabilities of not encountering (last column of m-array)
  for (k in 1:((nyears-1)*ns)){
    pi[k,nyears*ns-(ns-1)] <- 1-sum(pi[k,1:(nyears*ns-(ns-1))])
  } #t
  return(pi)
  returnType(double(2))
})



# Model with time-constant recovery parameters

# Write NIMBLE model file
Code <- nimbleCode({

# -------------------------------------------------
# Parameters:
# -------------------------------------------------
  # FEMALES
  # States (S):  
  # 1 Juvenile alive
  # 2 Not yet breeding at age 1 year
  # 3 Philopatric not yet breeding at age 2 year (inside the study area)
  # 4 Philopatric breeder (inside the study area)
  # 5 Philopatric breeder (inside the study area) that skips reproduction
  # 6 Non-philopatric not yet breeding at age 2 year (inside the study area)
  # 7 Non-philopatric breeder (inside the study area)
  # 8 Non-philopatric breeder (inside the study area) that skips reproduction
  # 9 Alive outside study area
  # 10 Recently dead by road-kill
  # 11 Recently dead by other causes
  # 12 Long-time dead (absorbing)
  # Observations (O):
  # 1 Seen as juvenile alive
  # 2 Seen as non-breeder at age 1 year
  # 3 Seen as philopatric non-breeder at age 2 year (inside the study area)
  # 4 Seen as philopatric breeder (inside the study area)
  # 5 Seen as philopatric non-breeder (inside the study area) (i.e. experienced breeder that skips reproduction)
  # 6 Seen as non-philopatric non-breeder at age 2 year (inside the study area)
  # 7 Seen as non-philopatric breeder (inside the study area)
  # 8 Seen as non-philopatric non-breeder (inside the study area) (i.e. experienced breeder that skips reproduction)
  # 9 Seen alive outside study area
  # 10 Recovered dead by road-kill
  # 11 Recovered dead by other causes
  # 12 Unobservable
  # MALES
  # States (S):
  # 1 Juvenile alive
  # 2 Alive at age 1 year
  # 3 Philopatric alive at age 2 year (inside the study area)
  # 4 Philopatric alive at age >=3 year (inside the study area)
  # 5 Non-philopatric at age 2 year (inside the study area)
  # 6 Non-philopatric at age >=3 year (inside the study area)
  # 7 Alive outside study area
  # 8 Recently dead by road-kill 
  # 9 Recently dead by other causes
  # 10 Long-time dead (absorbing)
  # Observations (O):
  # 1 Seen as juvenile alive
  # 2 Seen alive at age 1 year
  # 3 Seen as philopatric alive at age 2 year (inside the study area)
  # 4 Seen as philopatric alive at age >=3 year (inside the study area)
  # 5 Seen as non-philopatric alive at age 2 year (inside the study area)
  # 6 Seen as non-philopatric alive at age >=3 year (inside the study area)
  # 7 Seen alive outside study area
  # 8 Recovered dead by road-kill
  # 9 Recovered dead by other causes
  # 10 Unobservable
  # -------------------------------------------------
  # Priors and constraints
  # FIRST YEAR
  #============
  # FEMALES
  #~~~~~~~~~~
  # Mortality
  mr1F[1] <- mean.mr1F
  mo1F[1] <- mean.mo1F  
  log(mo2F[1]) <- mu.mo1F + ep.mo1F[1] # Background mortality hazard rate subadult female   
  log(mr2F[1]) <- mu.mr2F + beta*(NtH[1]-mean.C)/10 + ep.mr2F[1] # road-kill mortality hazard rate subadult female
  mr3F[1] <- mean.mr3F
  mo3F[1] <- mean.mo3F  
  # Survival
  s1F[1] <- exp(-(mr1F[1] + mo1F[1])) 
  s2F[1] <- exp(-(mr2F[1] + mo2F[1]))
  s3F[1] <- exp(-(mr3F[1] + mo3F[1]))
  # Roadkill mortality
  rk1F[1] <- (1 - s1F[1]) * (mr1F[1] / (mr1F[1] + mo1F[1]))
  rk2F[1] <- (1 - s2F[1]) * (mr2F[1] / (mr2F[1] + mo2F[1]))
  rk3F[1] <- (1 - s3F[1]) * (mr3F[1] / (mr3F[1] + mo3F[1]))
  # Background mortality
  o1F[1]  <- (1 - s1F[1]) * (mo1F[1] / (mr1F[1] + mo1F[1]))
  o2F[1]  <- (1 - s2F[1]) * (mo2F[1] / (mr2F[1] + mo2F[1]))
  o3F[1]  <- (1 - s3F[1]) * (mo3F[1] / (mr3F[1] + mo3F[1]))
  # Random effects
  ep.mr1F[1] <- 0
  ep.mo1F[1] <- 0	
  ep.mr2F[1] ~ T(dnorm(0, sd=sigma.mr2F),-10,10)
  ep.mo2F[1] ~ T(dnorm(0, sd=sigma.mo2F),-10,10)
  ep.mr3F[1] <- 0
  ep.mo3F[1] <- 0  
  # MALES
  #~~~~~~~~~~ 
  # Mortality
  mr1M[1] <- mean.mr1M
  mo1M[1] <- mean.mo1M 
  log(mr2M[1]) <- mu.mr2M + ep.mr2M[1] # road-kill mortality hazard rate subadult male
  log(mo2M[1]) <- mu.mo2M + ep.mo2M[1] # Background mortality hazard rate subadult male
  mr3M[1] <- mean.mr3M
  mo3M[1] <- mean.mo3M
  # Survival 
  s1M[1] <- exp(-(mr1M[1] + mo1M[1])) 
  s2M[1] <- exp(-(mr2M[1] + mo2M[1]))
  s3M[1] <- exp(-(mr3M[1] + mo3M[1]))
  # Roadkill mortality  
  rk1M[1] <- (1 - s1M[1]) * (mr1M[1] / (mr1M[1] + mo1M[1]))
  rk2M[1] <- (1 - s2M[1]) * (mr2M[1] / (mr2M[1] + mo2M[1]))
  rk3M[1] <- (1 - s3M[1]) * (mr3M[1] / (mr3M[1] + mo3M[1]))
  # Background mortality
  o1M[1]  <- (1 - s1M[1]) * (mo1M[1] / (mr1M[1] + mo1M[1]))
  o2M[1]  <- (1 - s2M[1]) * (mo2M[1] / (mr2M[1] + mo2M[1]))
  o3M[1]  <- (1 - s3M[1]) * (mo3M[1] / (mr3M[1] + mo3M[1]))
  # Random effects
  ep.mr1M[1] <- 0
  ep.mo1M[1] <- 0
  ep.mr2M[1] ~ T(dnorm(0, sd=sigma.mr2M),-10,10)
  ep.mo2M[1] ~ T(dnorm(0, sd=sigma.mo2M),-10,10)
  ep.mr3M[1] <- 0
  ep.mo3M[1] <- 0
  
  # SECOND AND FOLLOWING YEARS
  #============================
  for (t in 2:(n.occasions-1+K)){  
    # FEMALES
    #~~~~~~~~~~
    log(mr1F[t]) <- mu.mr1F + ep.mr1F[t] # road-kill mortality hazard rate juv
    log(mo1F[t]) <- mu.mo1F + ep.mo1F[t] # Background mortality hazard rate juv    
    log(mr2F[t]) <- mu.mr2F + beta*(NtH[t]-mean.C)/10 + ep.mr2F[t] # road-kill mortality hazard rate subadult
    log(mo2F[t]) <- mu.mo2F + ep.mo2F[t] # Background mortality hazard rate subadult
    log(mr3F[t]) <- mu.mr3F + ep.mr3F[t] # road-kill mortality hazard rate adult
    log(mo3F[t]) <- mu.mo3F + ep.mo3F[t] # Background mortality hazard rate adult
    # Survival	   
	s1F[t] <- exp(-(mr1F[t] + mo1F[t])) 
    s2F[t] <- exp(-(mr2F[t] + mo2F[t]))
    s3F[t] <- exp(-(mr3F[t] + mo3F[t]))
    # Roadkill mortality  	
    rk1F[t] <- (1 - s1F[t]) * (mr1F[t] / (mr1F[t] + mo1F[t]))
    rk2F[t] <- (1 - s2F[t]) * (mr2F[t] / (mr2F[t] + mo2F[t]))
    rk3F[t] <- (1 - s3F[t]) * (mr3F[t] / (mr3F[t] + mo3F[t]))
    # Background mortality    
	o1F[t]  <- (1 - s1F[t]) * (mo1F[t] / (mr1F[t] + mo1F[t]))
    o2F[t]  <- (1 - s2F[t]) * (mo2F[t] / (mr2F[t] + mo2F[t]))
    o3F[t]  <- (1 - s3F[t]) * (mo3F[t] / (mr3F[t] + mo3F[t]))
    # Random effects	
    ep.mr1F[t] ~ T(dnorm(0, sd=sigma.mr1F),-10,10)
    ep.mo1F[t] ~ T(dnorm(0, sd=sigma.mo1F),-10,10)	
    ep.mr2F[t] ~ T(dnorm(0, sd=sigma.mr2F),-10,10)
    ep.mo2F[t] ~ T(dnorm(0, sd=sigma.mo2F),-10,10)
    ep.mr3F[t] ~ T(dnorm(0, sd=sigma.mr3F),-10,10)
    ep.mo3F[t] ~ T(dnorm(0, sd=sigma.mo3F),-10,10)

    # MALES
    #~~~~~~~~~~
    # Survival
    log(mr1M[t]) <- mu.mr1M + ep.mr1M[t] # road-kill mortality hazard rate juv
    log(mo1M[t]) <- mu.mo1M + ep.mo1M[t] # Background mortality hazard rate juv   
    log(mr2M[t]) <- mu.mr2M + ep.mr2M[t] # road-kill mortality hazard rate subadult
    log(mo2M[t]) <- mu.mo2M + ep.mo2M[t] # Background mortality hazard rate subadult
    log(mr3M[t]) <- mu.mr3M + ep.mr3M[t] # road-kill mortality hazard rate adult
    log(mo3M[t]) <- mu.mo3M + ep.mo3M[t] # Background mortality hazard rate adult    
    # Survival	    
    s1M[t] <- exp(-(mr1M[t] + mo1M[t])) 
    s2M[t] <- exp(-(mr2M[t] + mo2M[t]))
    s3M[t] <- exp(-(mr3M[t] + mo3M[t]))
    # Roadkill mortality  
    rk1M[t] <- (1 - s1M[t]) * (mr1M[t] / (mr1M[t] + mo1M[t]))
    rk2M[t] <- (1 - s2M[t]) * (mr2M[t] / (mr2M[t] + mo2M[t]))
    rk3M[t] <- (1 - s3M[t]) * (mr3M[t] / (mr3M[t] + mo3M[t]))
    # Background mortality    
    o1M[t]  <- (1 - s1M[t]) * (mo1M[t] / (mr1M[t] + mo1M[t]))
    o2M[t]  <- (1 - s2M[t]) * (mo2M[t] / (mr2M[t] + mo2M[t]))
    o3M[t]  <- (1 - s3M[t]) * (mo3M[t] / (mr3M[t] + mo3M[t]))
    # Random effects	
    ep.mr1M[t] ~ T(dnorm(0, sd=sigma.mr1M),-10,10)
    ep.mo1M[t] ~ T(dnorm(0, sd=sigma.mo1M),-10,10)	
    ep.mr2M[t] ~ T(dnorm(0, sd=sigma.mr2M),-10,10)
    ep.mo2M[t] ~ T(dnorm(0, sd=sigma.mo2M),-10,10)
    ep.mr3M[t] ~ T(dnorm(0, sd=sigma.mr3M),-10,10)
    ep.mo3M[t] ~ T(dnorm(0, sd=sigma.mo3M),-10,10)
  }
  
  # Philopatry (or remains at release area) and fidelity
  for (t in 1:(n.occasions-1+K)){ 
    for(s in 1:2){
      # Non-philopatry/not remain in release area
      nonPhiF[s,t] <- mean.nonPhiF[s]  # Females that moved to an area other than the natal (or released) area      
      nonPhiM[s,t] <- mean.nonPhiM[s]  # Males that moved to an area other than the natal (or released) area
      # Fidelity
      FF[s,t] <- mean.FF[s]   # Female fidelity to the study area
      FM[s,t] <- mean.FM[s]   # Male fidelity to the study area
    }
  }
  
  for(s in 1:2){ # 1: wild-born; 2: released
    # Non-philopatry/not remain in release area
    mean.nonPhiF[s] ~ dunif(0, 1)   # Females. Not remain at release/born area
    mu.nonPhiF[s] <- log(mean.nonPhiF[s] / (1 - mean.nonPhiF[s]))
    mean.nonPhiM[s] ~ dunif(0, 1)   # Males. Not remain at release/born area
    mu.nonPhiM[s] <- log(mean.nonPhiM[s] / (1 - mean.nonPhiM[s]))		
    # Fidelity
    mean.FF[s] ~ dunif(0, 1)  # Females; 1: wild-born; 2: released
    mean.FM[s] ~ dunif(0, 1)  # Males; 1: wild-born; 2: released
  } 
  
  # Recovery and detection probabilities
  for(s in 1:2){ # 1: wild-born; 2: released
    for (t in 1:(n.occasions-1+K)){ 
      # Detection
      logit(p[s,t]) <- mu.p[s] + ep.p[s,t]    
      ep.p[s,t] ~ T(dnorm(0, sd=sigma.p[s]),-10,10)      
      # Recovery
      rF[s,t] <- mean.rF[s]
      rM[s,t] <- mean.rM[s]
    }	  
  }

  for(s in 1:2){ # 1: wild-born; 2: released
    # Detection
    mean.p[s] ~ dunif(0, 1) 
    mu.p[s] <- log(mean.p[s] /(1-mean.p[s]))
    sigma.p[s] ~ dunif(0,5)  
    
    # Recovery
    mean.rF[s] ~ dunif(0, 1)
    mean.rM[s] ~ dunif(0, 1)
  }  
  
  # Reproduction parameters
  ep.delta[1] <- 0   
  ep.alpha[1] <- 0
  ep.fec[1]   <- 0
  delta[1] <- mean.delta
  alpha[1] <- mean.alpha
  fec[1]   <- mean.fec
  for (t in 2:(n.occasions-1+K)){	
    logit(delta[t]) <- mu.delta + ep.delta[t]
    ep.delta[t] ~ T(dnorm(0, sd=sigma.delta),-10,10)
    alpha[t] <- mean.alpha
  }
  
  for (t in 2:n.occasions+K){		
    log(fec[t]) <- mu.fec + ep.fec[t]
    ep.fec[t]   ~  T(dnorm(0, sd=sigma.fec),-10,10)
  }
  
  beta ~ dunif(-5,5)
  
  mu.mr1F <- log(mean.mr1F)     
  mean.mr1F ~ dunif(0.001, 3)  # median road kill hazard rate
  sigma.mr1F ~ T(dnorm(0, 0.1),0,)

  mu.mr1M <- log(mean.mr1M)     
  mean.mr1M ~ dunif(0.001, 3)  # median road kill hazard rate
  sigma.mr1M ~ T(dnorm(0, 0.1),0,)
  
  mu.mo1F <- log(mean.mo1F)
  mean.mo1F ~ dunif(0.001, 3)
  sigma.mo1F ~ T(dnorm(0, 0.1),0,)

  mu.mo1M <- log(mean.mo1M)
  mean.mo1M ~ dunif(0.001, 3)
  sigma.mo1M ~ T(dnorm(0, 0.1),0,)
  
  mu.mr2F <- log(mean.mr2F)
  mean.mr2F ~ dunif(0.001, 3)
  sigma.mr2F ~ dunif(0,5)

  mu.mr2M <- log(mean.mr2M)
  mean.mr2M ~ dunif(0.001, 3)
  sigma.mr2M ~ dunif(0,5)
  
  mu.mo2F <- log(mean.mo2F)
  mean.mo2F ~ dunif(0.001, 3)
  sigma.mo2F ~ dunif(0,5)

  mu.mo2M <- log(mean.mo2M)
  mean.mo2M ~ dunif(0.001, 3)
  sigma.mo2M ~ dunif(0,5)   
  
  mu.mr3F <- log(mean.mr3F)
  mean.mr3F ~ dunif(0.001, 3)
  sigma.mr3F ~ dunif(0,5)

  mu.mr3M <- log(mean.mr3M)
  mean.mr3M ~ dunif(0.001, 3)
  sigma.mr3M ~ dunif(0,5)
  
  mu.mo3F <- log(mean.mo3F)
  mean.mo3F ~ dunif(0.001, 3)
  sigma.mo3F ~ dunif(0,5)

  mu.mo3M <- log(mean.mo3M)
  mean.mo3M ~ dunif(0.001, 3)
  sigma.mo3M ~ dunif(0,5)
  
  mean.alpha ~ dunif(0, 1)
  mu.alpha <- log(mean.alpha /(1-mean.alpha))
  # sigma.alpha ~ dunif(0,5)

  mean.delta ~ dunif(0, 1)
  mu.delta <- log(mean.delta / (1 - mean.delta))
  sigma.delta ~ dunif(0,5)
  
  mean.fec ~ dunif(0, 10)
  mu.fec <- log(mean.fec)
  sigma.fec ~ dunif(0,5)

  # Counts observation error
  sigma.y1 ~ dunif(0.5, 100)  # lower bound >0
  sigma.y2 ~ dunif(0.5, 100) 
 
  # Sex ratio of newborn individuals
  xi ~ dbeta(Jf, Jm)
  di ~ dbeta(R1, R2)
 
  # Initial population sizes
  n1wb[1] <- 0 
  n2wb[1] <- 0 
  n2rs[1] <- 0
  n2im[1] <- 0
  f2wb[1] <- 0
  f2rs[1] <- 0  
  f2im[1] <- 0  
  emwbF[1] <- 0
  emsF[1]  <- 0
  B3[1]   <- 0
  N.nonPhiwb[1] <- 0
  N.nonPhiRs[1] <- 0
  N.Phiwb[1] <- 0
  N.PhiRs[1] <- 0
  
  # Likelihood for population population count data (state-space model)
  # System process
  
  # FEMALES
  #=========
  for (t in 1:(n.occasions-1+K)){
    n1wb[t+1] ~ dpois(mu1F[t])                                  # n1wb: 1-year old wild-born females in year t
    mu1F[t] <- s1F[t]*fec[t]*xi*(F2[t] + delta[t]*B3[t])  
    
    # 2-YEARS OLD-NOT YET BREEDING
    n2wb[t+1] ~ dbin(s2F[t] * (1-alpha[t]) * FF[1,t], n1wb[t])  # 2-years old wild born females in year t that have never reproduced so far
    n2rs[t+1] ~ dbin(s2F[t] * (1-alpha[t]) * FF[2,t], RsF[t])   # 2-years old released females in year t that have never reproduced so far	
    n2im[t+1] ~ dbin(s2F[t] * (1-alpha[t]), round(IF[t]))       # 2-years old immigrant females in year t that have never reproduced so far 

    # 2-YEARS OLD-FIRST TIME BREEDING
    f2wb[t+1] ~ dbin(s2F[t] * alpha[t] * FF[1,t], n1wb[t])      # 2-years old wild born females in year t that reproduced for the first time in t
    f2rs[t+1] ~ dbin(s2F[t] * alpha[t] * FF[2,t], RsF[t])       # 2-years old released females in year t that reproduced for the first time in t
    f2im[t+1] ~ dbin(s2F[t] * alpha[t], round(IF[t]))           # 2-years old immigrant females in year t that reproduced for the first time in t

    # EMIGRANTS 
    emwbF[t+1] ~ dbin(s2F[t]*(1-FF[1,t]), n1wb[t])              # 2-years old wild born emigrant females in year t
    emsF[t+1] ~ dbin(s2F[t]*(1-FF[2,t]), RsF[t])                # 2-years old emigrant released in year t

    # BREEDERS
    B3[t+1] ~ dbin(s3F[t], N2T[t] + B3[t])                      # 3-year or older in year t
	
    # 2-YEARS NON-PHILOPATRIC (WILD BORN) - NOT REMAIN AT RELEASE AREA (RELEASED)
    N.nonPhiwb[t+1] ~ dbin(s2F[t] * FF[1,t] * nonPhiF[1,t], n1wb[t])     # wild-born
    N.nonPhiRs[t+1] ~ dbin(s2F[t] * FF[2,t] * nonPhiF[2,t], RsF[t])	     # release
	
    # 2-YEARS PHILOPATRIC (WILD BORN)- REMAIN AT RELEASE AREA (RELEASED)	
    N.Phiwb[t+1] ~ dbin(s2F[t] * FF[1,t] * (1-nonPhiF[1,t]), n1wb[t])    # wild-born
    N.PhiRs[t+1] ~ dbin(s2F[t] * FF[2,t]* (1-nonPhiF[2,t]), RsF[t])	     # release
  }  
  
  # IMMIGRATION FEMALES
  IF[1] <- 0
  log.IF[1] <- 0
  for(t in 2:n.occasions){
    log.IF[t] ~ dnorm(0, sd=10)
    IF[t] <- exp(log.IF[t])
	CimF[t] ~ dpois(IF[t])
  }
  mean.IF <- mean(IF[2:n.occasions])
  sigma.IF <- sd(IF[2:n.occasions])             # we assumed the same immigration
  for (t in (n.occasions+1):(n.occasions+K)) { 
    IF[t] ~ T(dnorm(mean.IF, sd=sigma.IF),0,)
  }
  
  
  # SSM FEMALES
  FRep[1] <- F2[1] + B3[1]*mean.delta 
  for (t in 2:n.occasions){
    FRep[t] <- F2[t] + B3[t]*delta[t-1]
    # Observation process
    C.FRep[t] ~ dnorm(FRep[t], sd=sigma.y1)	
  }
  for(t in (n.occasions+1):(n.occasions+K)){
    FRep[t] <- F2[t] + B3[t]*delta[t-1]
  }

  # Derived quantities
  for (t in 1:(n.occasions+K)) {
    N1[t] <- n1wb[t] + RsF[t] +IF[t]        # 1-year old
    N2[t] <- n2wb[t] + n2rs[t] + n2im[t]    # 2-years old that have never reproduced
    F2[t] <- f2wb[t] + f2rs[t] + f2im[t]    # 2-years old that reproduced for the first time
    N2T[t] <- N2[t]+F2[t]
    EMF[t] <- emwbF[t] + emsF[t]            # Emigrants  
    # Total females
    NtH[t]<- N1[t] + N2T[t] + B3[t]	  
  } 
  
  # MALES
  #=========
  n1wbM[1] <- 0
  N3M[1] <- 0
  N.nonPhiwbM[1] <- 0
  N.nonPhiRsM[1] <- 0
  N.PhiwbM[1] <- 0
  N.PhiRsM[1] <- 0
  emwbM[1] <- 0 # 2-years old wild born emigrants in year 1
  emsM[1]  <- 0 # 2-years old released emigrants in year 1
  nim[1] <- 0
  
  for (t in 1:(n.occasions-1+K)){
    n1wbM[t+1] ~ dpois(mu1M[t])                                         # n1wbM: 1-year old wild-born males in year t
    mu1M[t] <- s1M[t]*fec[t]*(1-xi)*(F2[t] + delta[t]*B3[t])
    
    # 2-YEARS OLD
    # 2-years non-philopatric (wild born) / not remain at release area (released)
    N.nonPhiwbM[t+1] ~ dbin(s2M[t] * FM[1,t] * nonPhiM[1,t], n1wbM[t])  # wild-born
    N.nonPhiRsM[t+1] ~ dbin(s2M[t] * FM[2,t] * nonPhiM[2,t], RsM[t])    # release
	
    # 2-years philopatric (wild born) / remain at release area (released)	
    N.PhiwbM[t+1] ~ dbin(s2M[t] * FM[1,t] * (1-nonPhiM[1,t]), n1wbM[t]) # wild-born
    N.PhiRsM[t+1] ~ dbin(s2M[t] * FM[2,t]* (1-nonPhiM[2,t]), RsM[t])	# release

    # Immigrants
    nim[t+1] ~ dbin(s2M[t], round(IM[t]))

    # EMIGRANTS 
    emwbM[t+1] ~ dbin(s2M[t]*(1-FM[1,t]), n1wbM[t])       # 2-years old wild born emigrants in year t
    emsM[t+1] ~ dbin(s2M[t]*(1-FM[2,t]), RsM[t])          # 2-years old released emigrants in year t

    # 3-YEARS AND OLDER
    N3M[t+1] ~ dbin(s3M[t], N2M[t] + N3M[t])
  }	
  
  # IMMIGRATION MALES
  IM[1]<-0
  log.IM[1] <- 0
  for(t in 2:n.occasions){
    log.IM[t] ~ dnorm(0, sd=10)
    IM[t] <- exp(log.IM[t])	
    CimM[t] ~ dpois(IM[t])
  }
  mean.IM <- mean(IM[2:n.occasions])
  sigma.IM <- sd(IM[2:n.occasions])             # we assumed the same immigration
  for (t in (n.occasions+1):(n.occasions+K)) { 
    IM[t] ~ T(dnorm(mean.IM, sd=sigma.IM),0,)
  }
 
  # SSM MALES
  for(t in 1:n.occasions){
    NMrep[t] <- N2M[t] + N3M[t]
    C.Mrep [t] ~ dnorm(NMrep[t], sd=sigma.y2)
  }

  # Derived quantities
  for (t in 1:(n.occasions+K)) {
    N1M[t] <- n1wbM[t] + RsM[t] + IM[t]          # 1-year old
    N2wb[t] <- N.nonPhiwbM[t] + N.PhiwbM[t]
    N2Rs[t] <- N.nonPhiRsM[t] + N.PhiRsM[t]
	N2M[t] <- N2wb[t] + N2Rs[t] + nim[t]     # 2-years old
    NtM[t]<- N1M[t] + N2M[t] + N3M[t]
    # Population size (F&M)
    NT[t] <- NtH[t] + NtM[t]
  }   
  
  # Likelihood for reproductive data: Poisson regression
  for (t in 1:n.occasions) {
    J[t] ~ dpois(rho[t])                   # Number of kittens
    log(rho[t]) <- log(R[t])+log(fec[t])   # R: Number of females that have raised these kittens
  }
  
  # Check whether the population is extinct in the future
  for (t in 1:K){
    extinct[t] <- equals(NT[n.occasions+t], 0)
  }
  
  # FEMALES
  #===========
  # Define state-transition and re-encounter probabilities
  for (t in 1:(n.occasions-1)){ 
      
      # Group 1 (wild born)
      ps1F[1,1,t] <- 0
      ps1F[1,2,t] <- s1F[t]
      ps1F[1,3,t] <- 0
      ps1F[1,4,t] <- 0
      ps1F[1,5,t] <- 0
      ps1F[1,6,t] <- 0
      ps1F[1,7,t] <- 0
      ps1F[1,8,t] <- 0
      ps1F[1,9,t] <- 0  
      ps1F[1,10,t] <- rk1F[t]
      ps1F[1,11,t] <- o1F[t]
      ps1F[1,12,t] <- 0
      
      ps1F[2,1,t] <- 0
      ps1F[2,2,t] <- 0
      ps1F[2,3,t] <- s2F[t] * FF[1,t] * (1-alpha[t]) * (1-nonPhiF[1,t])
      ps1F[2,4,t] <- s2F[t] * FF[1,t] * alpha[t] * (1-nonPhiF[1,t])
      ps1F[2,5,t] <- 0
      ps1F[2,6,t] <- s2F[t] * FF[1,t] * (1-alpha[t]) * nonPhiF[1,t]
      ps1F[2,7,t] <- s2F[t] * FF[1,t] * alpha[t] * nonPhiF[1,t]
      ps1F[2,8,t] <- 0
      ps1F[2,9,t] <- s2F[t] * (1-FF[1,t])
      ps1F[2,10,t] <- rk2F[t] 
      ps1F[2,11,t] <- o2F[t]
      ps1F[2,12,t] <- 0
      
      ps1F[3,1,t] <- 0
      ps1F[3,2,t] <- 0
      ps1F[3,3,t] <- 0
      ps1F[3,4,t] <- s3F[t] * delta[t]
      ps1F[3,5,t] <- s3F[t] * (1-delta[t])
      ps1F[3,6,t] <- 0
      ps1F[3,7,t] <- 0
      ps1F[3,8,t] <- 0
      ps1F[3,9,t] <- 0
      ps1F[3,10,t] <- rk3F[t]
      ps1F[3,11,t] <- o3F[t] 
      ps1F[3,12,t] <- 0
      
      ps1F[4,1,t] <- 0
      ps1F[4,2,t] <- 0
      ps1F[4,3,t] <- 0 
      ps1F[4,4,t] <- s3F[t] * delta[t]
      ps1F[4,5,t] <- s3F[t] * (1-delta[t])
      ps1F[4,6,t] <- 0
      ps1F[4,7,t] <- 0
      ps1F[4,8,t] <- 0   
      ps1F[4,9,t] <- 0
      ps1F[4,10,t] <- rk3F[t]
      ps1F[4,11,t] <- o3F[t]
      ps1F[4,12,t] <- 0
    
      ps1F[5,1,t] <- 0
      ps1F[5,2,t] <- 0
      ps1F[5,3,t] <- 0 
      ps1F[5,4,t] <- s3F[t] * delta[t]
      ps1F[5,5,t] <- s3F[t] * (1-delta[t])
      ps1F[5,6,t] <- 0
      ps1F[5,7,t] <- 0
      ps1F[5,8,t] <- 0 
      ps1F[5,9,t] <- 0
      ps1F[5,10,t] <- rk3F[t]
      ps1F[5,11,t] <- o3F[t]
      ps1F[5,12,t] <- 0
      
      ps1F[6,1,t] <- 0
      ps1F[6,2,t] <- 0
      ps1F[6,3,t] <- 0
      ps1F[6,4,t] <- 0
      ps1F[6,5,t] <- 0
      ps1F[6,6,t] <- 0
      ps1F[6,7,t] <- s3F[t] * delta[t]
      ps1F[6,8,t] <- s3F[t] * (1-delta[t])
      ps1F[6,9,t] <- 0
      ps1F[6,10,t] <- rk3F[t]
      ps1F[6,11,t] <- o3F[t]
      ps1F[6,12,t] <- 0
      
      ps1F[7,1,t] <- 0
      ps1F[7,2,t] <- 0
      ps1F[7,3,t] <- 0
      ps1F[7,4,t] <- 0
      ps1F[7,5,t] <- 0
      ps1F[7,6,t] <- 0
      ps1F[7,7,t] <- s3F[t] * delta[t]
      ps1F[7,8,t] <- s3F[t] * (1-delta[t])
      ps1F[7,9,t] <- 0
      ps1F[7,10,t] <- rk3F[t]
      ps1F[7,11,t] <- o3F[t]
      ps1F[7,12,t] <- 0
    
      ps1F[8,1,t] <- 0
      ps1F[8,2,t] <- 0
      ps1F[8,3,t] <- 0
      ps1F[8,4,t] <- 0
      ps1F[8,5,t] <- 0
      ps1F[8,6,t] <- 0
      ps1F[8,7,t] <- s3F[t] * delta[t]
      ps1F[8,8,t] <- s3F[t] * (1-delta[t])
      ps1F[8,9,t] <- 0
      ps1F[8,10,t] <- rk3F[t]
      ps1F[8,11,t] <- o3F[t]
      ps1F[8,12,t] <- 0
    
      ps1F[9,1,t] <- 0
      ps1F[9,2,t] <- 0
      ps1F[9,3,t] <- 0
      ps1F[9,4,t] <- 0
      ps1F[9,5,t] <- 0
      ps1F[9,6,t] <- 0
      ps1F[9,7,t] <- 0
      ps1F[9,8,t] <- 0
      ps1F[9,9,t] <- s3F[t]
      ps1F[9,10,t] <- rk3F[t]
      ps1F[9,11,t] <- o3F[t] 
      ps1F[9,12,t] <- 0
      
      ps1F[10,1,t] <- 0
      ps1F[10,2,t] <- 0
      ps1F[10,3,t] <- 0
      ps1F[10,4,t] <- 0
      ps1F[10,5,t] <- 0
      ps1F[10,6,t] <- 0
      ps1F[10,7,t] <- 0
      ps1F[10,8,t] <- 0
      ps1F[10,9,t] <- 0
      ps1F[10,10,t] <- 0
      ps1F[10,11,t] <- 0
      ps1F[10,12,t] <- 1
    
      ps1F[11,1,t] <- 0
      ps1F[11,2,t] <- 0
      ps1F[11,3,t] <- 0
      ps1F[11,4,t] <- 0
      ps1F[11,5,t] <- 0
      ps1F[11,6,t] <- 0
      ps1F[11,7,t] <- 0
      ps1F[11,8,t] <- 0
      ps1F[11,9,t] <- 0
      ps1F[11,10,t] <- 0
      ps1F[11,11,t] <- 0
      ps1F[11,12,t] <- 1
      
      ps1F[12,1,t] <- 0
      ps1F[12,2,t] <- 0
      ps1F[12,3,t] <- 0
      ps1F[12,4,t] <- 0
      ps1F[12,5,t] <- 0
      ps1F[12,6,t] <- 0
      ps1F[12,7,t] <- 0
      ps1F[12,8,t] <- 0
      ps1F[12,9,t] <- 0
      ps1F[12,10,t] <- 0
      ps1F[12,11,t] <- 0  
      ps1F[12,12,t] <- 1
    
      # Group 2 (release)
      ps2F[1,1,t] <- 0
      ps2F[1,2,t] <- s1F[t]
      ps2F[1,3,t] <- 0
      ps2F[1,4,t] <- 0
      ps2F[1,5,t] <- 0
      ps2F[1,6,t] <- 0
      ps2F[1,7,t] <- 0
      ps2F[1,8,t] <- 0
      ps2F[1,9,t] <- 0  
      ps2F[1,10,t] <- rk1F[t]
      ps2F[1,11,t] <- o1F[t]
      ps2F[1,12,t] <- 0
      
      ps2F[2,1,t] <- 0
      ps2F[2,2,t] <- 0
      ps2F[2,3,t] <- s2F[t] * FF[2,t] * (1-alpha[t]) * (1-nonPhiF[2,t])
      ps2F[2,4,t] <- s2F[t] * FF[2,t] * alpha[t] * (1-nonPhiF[2,t])
      ps2F[2,5,t] <- 0
      ps2F[2,6,t] <- s2F[t] * FF[2,t] * (1-alpha[t]) * nonPhiF[2,t]
      ps2F[2,7,t] <- s2F[t] * FF[2,t] * alpha[t] * nonPhiF[2,t]
      ps2F[2,8,t] <- 0
      ps2F[2,9,t] <- s2F[t] * (1-FF[2,t])
      ps2F[2,10,t] <- rk2F[t] 
      ps2F[2,11,t] <- o2F[t] 
      ps2F[2,12,t] <- 0
      
      ps2F[3,1,t] <- 0
      ps2F[3,2,t] <- 0
      ps2F[3,3,t] <- 0
      ps2F[3,4,t] <- s3F[t] * delta[t]
      ps2F[3,5,t] <- s3F[t] * (1-delta[t])
      ps2F[3,6,t] <- 0
      ps2F[3,7,t] <- 0
      ps2F[3,8,t] <- 0
      ps2F[3,9,t] <- 0
      ps2F[3,10,t] <- rk3F[t] 
      ps2F[3,11,t] <- o3F[t]  
      ps2F[3,12,t] <- 0
      
      ps2F[4,1,t] <- 0
      ps2F[4,2,t] <- 0
      ps2F[4,3,t] <- 0 
      ps2F[4,4,t] <- s3F[t] * delta[t]
      ps2F[4,5,t] <- s3F[t] * (1-delta[t])
      ps2F[4,6,t] <- 0
      ps2F[4,7,t] <- 0
      ps2F[4,8,t] <- 0   
      ps2F[4,9,t] <- 0
      ps2F[4,10,t] <- rk3F[t]
      ps2F[4,11,t] <- o3F[t]
      ps2F[4,12,t] <- 0
    
      ps2F[5,1,t] <- 0
      ps2F[5,2,t] <- 0
      ps2F[5,3,t] <- 0 
      ps2F[5,4,t] <- s3F[t] * delta[t]
      ps2F[5,5,t] <- s3F[t] * (1-delta[t])
      ps2F[5,6,t] <- 0
      ps2F[5,7,t] <- 0
      ps2F[5,8,t] <- 0 
      ps2F[5,9,t] <- 0
      ps2F[5,10,t] <- rk3F[t]
      ps2F[5,11,t] <- o3F[t]
      ps2F[5,12,t] <- 0
      
      ps2F[6,1,t] <- 0
      ps2F[6,2,t] <- 0
      ps2F[6,3,t] <- 0
      ps2F[6,4,t] <- 0
      ps2F[6,5,t] <- 0 
      ps2F[6,6,t] <- 0
      ps2F[6,7,t] <- s3F[t] * delta[t]
      ps2F[6,8,t] <- s3F[t] * (1-delta[t])
      ps2F[6,9,t] <- 0
      ps2F[6,10,t] <- rk3F[t]
      ps2F[6,11,t] <- o3F[t]
      ps2F[6,12,t] <- 0
      
      ps2F[7,1,t] <- 0
      ps2F[7,2,t] <- 0
      ps2F[7,3,t] <- 0
      ps2F[7,4,t] <- 0
      ps2F[7,5,t] <- 0
      ps2F[7,6,t] <- 0
      ps2F[7,7,t] <- s3F[t] * delta[t]
      ps2F[7,8,t] <- s3F[t] * (1-delta[t])
      ps2F[7,9,t] <- 0
      ps2F[7,10,t] <- rk3F[t]
      ps2F[7,11,t] <- o3F[t]
      ps2F[7,12,t] <- 0
    
      ps2F[8,1,t] <- 0
      ps2F[8,2,t] <- 0
      ps2F[8,3,t] <- 0
      ps2F[8,4,t] <- 0
      ps2F[8,5,t] <- 0
      ps2F[8,6,t] <- 0
      ps2F[8,7,t] <- s3F[t] * delta[t]
      ps2F[8,8,t] <- s3F[t] * (1-delta[t])
      ps2F[8,9,t] <- 0
      ps2F[8,10,t] <- rk3F[t]
      ps2F[8,11,t] <- o3F[t]
      ps2F[8,12,t] <- 0
    
      ps2F[9,1,t] <- 0
      ps2F[9,2,t] <- 0
      ps2F[9,3,t] <- 0
      ps2F[9,4,t] <- 0
      ps2F[9,5,t] <- 0
      ps2F[9,6,t] <- 0
      ps2F[9,7,t] <- 0
      ps2F[9,8,t] <- 0
      ps2F[9,9,t] <- s3F[t]
      ps2F[9,10,t] <- rk3F[t]
      ps2F[9,11,t] <- o3F[t] 
      ps2F[9,12,t] <- 0
      
      ps2F[10,1,t] <- 0
      ps2F[10,2,t] <- 0
      ps2F[10,3,t] <- 0
      ps2F[10,4,t] <- 0
      ps2F[10,5,t] <- 0
      ps2F[10,6,t] <- 0
      ps2F[10,7,t] <- 0
      ps2F[10,8,t] <- 0
      ps2F[10,9,t] <- 0
      ps2F[10,10,t] <- 0
      ps2F[10,11,t] <- 0
      ps2F[10,12,t] <- 1
    
      ps2F[11,1,t] <- 0
      ps2F[11,2,t] <- 0
      ps2F[11,3,t] <- 0
      ps2F[11,4,t] <- 0
      ps2F[11,5,t] <- 0
      ps2F[11,6,t] <- 0
      ps2F[11,7,t] <- 0
      ps2F[11,8,t] <- 0
      ps2F[11,9,t] <- 0
      ps2F[11,10,t] <- 0
      ps2F[11,11,t] <- 0
      ps2F[11,12,t] <- 1
      
      ps2F[12,1,t] <- 0
      ps2F[12,2,t] <- 0
      ps2F[12,3,t] <- 0
      ps2F[12,4,t] <- 0
      ps2F[12,5,t] <- 0
      ps2F[12,6,t] <- 0
      ps2F[12,7,t] <- 0
      ps2F[12,8,t] <- 0
      ps2F[12,9,t] <- 0
      ps2F[12,10,t] <- 0
      ps2F[12,11,t] <- 0  
      ps2F[12,12,t] <- 1
   

      # Define probabilities of O(t) given S(t)
      # Group 1 (wild-born)
      po1F[1,t] <- p[1,t]
      po1F[2,t] <- p[1,t]
      po1F[3,t] <- p[1,t]
      po1F[4,t] <- p[1,t]
      po1F[5,t] <- p[1,t]
      po1F[6,t] <- p[1,t]
      po1F[7,t] <- p[1,t]
      po1F[8,t] <- p[1,t]
      po1F[9,t] <- p[1,t]
      po1F[10,t] <- rF[1,t]
      po1F[11,t] <- rF[1,t]  
      po1F[12,t] <- 0
      
      # Group 2 (release)
      po2F[1,t] <- p[2,t]
      po2F[2,t] <- p[2,t]
      po2F[3,t] <- p[2,t]
      po2F[4,t] <- p[2,t]
      po2F[5,t] <- p[2,t]
      po2F[6,t] <- p[2,t]
      po2F[7,t] <- p[2,t]
      po2F[8,t] <- p[2,t]
      po2F[9,t] <- p[2,t]
      po2F[10,t] <- rF[2,t]
      po2F[11,t] <- rF[2,t]  
      po2F[12,t] <- 0


  }

  # MALES
  #=========== 
  # Define state-transition and re-encounter probabilities for males
  for (t in 1:(n.occasions-1)){ 
  
      # Group 1 (wild-born)
      ps1M[1,1,t] <- 0
      ps1M[1,2,t] <- s1M[t]
      ps1M[1,3,t] <- 0
      ps1M[1,4,t] <- 0
      ps1M[1,5,t] <- 0
      ps1M[1,6,t] <- 0 
      ps1M[1,7,t] <- 0 	  
      ps1M[1,8,t] <- rk1M[t]
      ps1M[1,9,t] <- o1M[t]
      ps1M[1,10,t] <- 0
      
      ps1M[2,1,t] <- 0
      ps1M[2,2,t] <- 0
      ps1M[2,3,t] <- s2M[t] * FM[1,t] * (1-nonPhiM[1,t])
      ps1M[2,4,t] <- 0
      ps1M[2,5,t] <- s2M[t] * FM[1,t] * nonPhiM[1,t]
      ps1M[2,6,t] <- 0	  
      ps1M[2,7,t] <- s2M[t] * (1-FM[1,t])
      ps1M[2,8,t] <- rk2M[t] 
      ps1M[2,9,t] <- o2M[t]
      ps1M[2,10,t] <- 0
      
      ps1M[3,1,t] <- 0
      ps1M[3,2,t] <- 0
      ps1M[3,3,t] <- 0
      ps1M[3,4,t] <- s3M[t]
      ps1M[3,5,t] <- 0
      ps1M[3,6,t] <- 0
      ps1M[3,7,t] <- 0
      ps1M[3,8,t] <- rk3M[t]
      ps1M[3,9,t] <- o3M[t] 
      ps1M[3,10,t] <- 0
      
      ps1M[4,1,t] <- 0
      ps1M[4,2,t] <- 0
      ps1M[4,3,t] <- 0 
      ps1M[4,4,t] <- s3M[t]
      ps1M[4,5,t] <- 0
      ps1M[4,6,t] <- 0
      ps1M[4,7,t] <- 0
      ps1M[4,8,t] <- rk3M[t]
      ps1M[4,9,t] <- o3M[t]
      ps1M[4,10,t] <- 0
    
      ps1M[5,1,t] <- 0
      ps1M[5,2,t] <- 0
      ps1M[5,3,t] <- 0 
      ps1M[5,4,t] <- 0
      ps1M[5,5,t] <- 0
      ps1M[5,6,t] <- s3M[t]
      ps1M[5,7,t] <- 0
      ps1M[5,8,t] <- rk3M[t]
      ps1M[5,9,t] <- o3M[t]
      ps1M[5,10,t] <- 0
     
      ps1M[6,1,t] <- 0
      ps1M[6,2,t] <- 0
      ps1M[6,3,t] <- 0 
      ps1M[6,4,t] <- 0
      ps1M[6,5,t] <- 0
      ps1M[6,6,t] <- s3M[t]
      ps1M[6,7,t] <- 0
      ps1M[6,8,t] <- rk3M[t]
      ps1M[6,9,t] <- o3M[t]
      ps1M[6,10,t] <- 0
    
      ps1M[7,1,t] <- 0
      ps1M[7,2,t] <- 0
      ps1M[7,3,t] <- 0 
      ps1M[7,4,t] <- 0
      ps1M[7,5,t] <- 0
      ps1M[7,6,t] <- 0
      ps1M[7,7,t] <- s3M[t]
      ps1M[7,8,t] <- rk3M[t]
      ps1M[7,9,t] <- o3M[t]
      ps1M[7,10,t] <- 0
	  
      ps1M[8,1,t] <- 0
      ps1M[8,2,t] <- 0
      ps1M[8,3,t] <- 0
      ps1M[8,4,t] <- 0
      ps1M[8,5,t] <- 0
      ps1M[8,6,t] <- 0
      ps1M[8,7,t] <- 0	  
      ps1M[8,8,t] <- 0
      ps1M[8,9,t] <- 0
      ps1M[8,10,t] <- 1
    
      ps1M[9,1,t] <- 0
      ps1M[9,2,t] <- 0
      ps1M[9,3,t] <- 0
      ps1M[9,4,t] <- 0
      ps1M[9,5,t] <- 0
      ps1M[9,6,t] <- 0
      ps1M[9,7,t] <- 0
      ps1M[9,8,t] <- 0
      ps1M[9,9,t] <- 0
      ps1M[9,10,t] <- 1
      
      ps1M[10,1,t] <- 0
      ps1M[10,2,t] <- 0
      ps1M[10,3,t] <- 0
      ps1M[10,4,t] <- 0
      ps1M[10,5,t] <- 0
      ps1M[10,6,t] <- 0
      ps1M[10,7,t] <- 0 
      ps1M[10,8,t] <- 0
      ps1M[10,9,t] <- 0
      ps1M[10,10,t] <- 1
    
      # Group 2 (release)
      ps2M[1,1,t] <- 0
      ps2M[1,2,t] <- s1M[t]
      ps2M[1,3,t] <- 0
      ps2M[1,4,t] <- 0
      ps2M[1,5,t] <- 0
      ps2M[1,6,t] <- 0 
      ps2M[1,7,t] <- 0 	  
      ps2M[1,8,t] <- rk1M[t]
      ps2M[1,9,t] <- o1M[t]
      ps2M[1,10,t] <- 0
      
      ps2M[2,1,t] <- 0
      ps2M[2,2,t] <- 0
      ps2M[2,3,t] <- s2M[t] * FM[2,t] * (1-nonPhiM[2,t])
      ps2M[2,4,t] <- 0
      ps2M[2,5,t] <- s2M[t] * FM[2,t] * nonPhiM[2,t]
      ps2M[2,6,t] <- 0
      ps2M[2,7,t] <- s2M[t] * (1-FM[2,t])
      ps2M[2,8,t] <- rk2M[t] 
      ps2M[2,9,t] <- o2M[t]
      ps2M[2,10,t] <- 0
      
      ps2M[3,1,t] <- 0
      ps2M[3,2,t] <- 0
      ps2M[3,3,t] <- 0
      ps2M[3,4,t] <- s3M[t]
      ps2M[3,5,t] <- 0
      ps2M[3,6,t] <- 0
      ps2M[3,7,t] <- 0
      ps2M[3,8,t] <- rk3M[t]
      ps2M[3,9,t] <- o3M[t] 
      ps2M[3,10,t] <- 0
      
      ps2M[4,1,t] <- 0
      ps2M[4,2,t] <- 0
      ps2M[4,3,t] <- 0 
      ps2M[4,4,t] <- s3M[t]
      ps2M[4,5,t] <- 0
      ps2M[4,6,t] <- 0
      ps2M[4,7,t] <- 0
      ps2M[4,8,t] <- rk3M[t]
      ps2M[4,9,t] <- o3M[t]
      ps2M[4,10,t] <- 0
    
      ps2M[5,1,t] <- 0
      ps2M[5,2,t] <- 0
      ps2M[5,3,t] <- 0 
      ps2M[5,4,t] <- 0
      ps2M[5,5,t] <- 0
      ps2M[5,6,t] <- s3M[t]
      ps2M[5,7,t] <- 0
      ps2M[5,8,t] <- rk3M[t]
      ps2M[5,9,t] <- o3M[t]
      ps2M[5,10,t] <- 0
     
      ps2M[6,1,t] <- 0
      ps2M[6,2,t] <- 0
      ps2M[6,3,t] <- 0 
      ps2M[6,4,t] <- 0
      ps2M[6,5,t] <- 0
      ps2M[6,6,t] <- s3M[t]
      ps2M[6,7,t] <- 0
      ps2M[6,8,t] <- rk3M[t]
      ps2M[6,9,t] <- o3M[t]
      ps2M[6,10,t] <- 0
    
      ps2M[7,1,t] <- 0
      ps2M[7,2,t] <- 0
      ps2M[7,3,t] <- 0 
      ps2M[7,4,t] <- 0
      ps2M[7,5,t] <- 0
      ps2M[7,6,t] <- 0
      ps2M[7,7,t] <- s3M[t]
      ps2M[7,8,t] <- rk3M[t]
      ps2M[7,9,t] <- o3M[t]
      ps2M[7,10,t] <- 0
	  
      ps2M[8,1,t] <- 0
      ps2M[8,2,t] <- 0
      ps2M[8,3,t] <- 0
      ps2M[8,4,t] <- 0
      ps2M[8,5,t] <- 0
      ps2M[8,6,t] <- 0
      ps2M[8,7,t] <- 0	  
      ps2M[8,8,t] <- 0
      ps2M[8,9,t] <- 0
      ps2M[8,10,t] <- 1
    
      ps2M[9,1,t] <- 0
      ps2M[9,2,t] <- 0
      ps2M[9,3,t] <- 0
      ps2M[9,4,t] <- 0
      ps2M[9,5,t] <- 0
      ps2M[9,6,t] <- 0
      ps2M[9,7,t] <- 0
      ps2M[9,8,t] <- 0
      ps2M[9,9,t] <- 0
      ps2M[9,10,t] <- 1
      
      ps2M[10,1,t] <- 0
      ps2M[10,2,t] <- 0
      ps2M[10,3,t] <- 0
      ps2M[10,4,t] <- 0
      ps2M[10,5,t] <- 0
      ps2M[10,6,t] <- 0
      ps2M[10,7,t] <- 0 
      ps2M[10,8,t] <- 0
      ps2M[10,9,t] <- 0
      ps2M[10,10,t] <- 1

      # Define probabilities of O(t) given S(t)
      # Group 1 (wild-born)
      po1M[1,t] <- p[1,t]
      po1M[2,t] <- p[1,t]
      po1M[3,t] <- p[1,t]
      po1M[4,t] <- p[1,t]
      po1M[5,t] <- p[1,t]
      po1M[6,t] <- p[1,t]
      po1M[7,t] <- p[1,t]
      po1M[8,t] <- rM[1,t]
      po1M[9,t] <- rM[1,t]  
      po1M[10,t] <- 0
      
      # Group 2 (release)
      po2M[1,t] <- p[2,t]
      po2M[2,t] <- p[2,t]
      po2M[3,t] <- p[2,t]
      po2M[4,t] <- p[2,t]
      po2M[5,t] <- p[2,t]
      po2M[6,t] <- p[2,t]
      po2M[7,t] <- p[2,t]
      po2M[8,t] <- rM[2,t]
      po2M[9,t] <- rM[2,t]  
      po2M[10,t] <- 0

  } 

  # FEMALES
  #~~~~~~~~~~~~
  # Calculate probability of non-encounter (dq) and reshape the array for the encounter
  # probabilities
  for (t in 1:(n.occasions-1)){
    for (s in 1:ns1){
      dp1F[s,s,t] <- po1F[s,t]
      dq1F[s,s,t] <- 1-po1F[s,t]
      dp2F[s,s,t] <- po2F[s,t]
      dq2F[s,s,t] <- 1-po2F[s,t]
    } #s
    for (s in 1:(ns1-1)){
      for (m in (s+1):ns1){
        dp1F[s,m,t] <- 0
        dq1F[s,m,t] <- 0
        dp2F[s,m,t] <- 0
        dq2F[s,m,t] <- 0
      } #s
    } #m
    for (s in 2:ns1){
      for (m in 1:(s-1)){
        dp1F[s,m,t] <- 0
        dq1F[s,m,t] <- 0
        dp2F[s,m,t] <- 0
        dq2F[s,m,t] <- 0
      } #s
    } #m
  } #t
  
  # Define the multinomial likelihood
  for (t in 1:((n.occasions-1)*ns1)){
    marrF[t,1:(n.occasions *ns1-(ns1-1)),1] ~ dmulti(pi1F[t,1:(n.occasions*ns1-(ns1-1))], released1F[t])
    marrF[t,1:(n.occasions *ns1-(ns1-1)),2] ~ dmulti(pi2F[t,1:(n.occasions*ns1-(ns1-1))], released2F[t])
  }

  # Define cell probabilities of the multistate m-array
  pi1F[1:((n.occasions-1)*ns1), 1:(n.occasions*ns1-(ns1-1))] <- marrProb(ps1F[1:ns1,1:ns1,1:(n.occasions-1)], dp1F[1:ns1,1:ns1,1:(n.occasions-1)], dq1F[1:ns1,1:ns1,1:(n.occasions-1)])
  pi2F[1:((n.occasions-1)*ns1), 1:(n.occasions*ns1-(ns1-1))] <- marrProb(ps2F[1:ns1,1:ns1,1:(n.occasions-1)], dp2F[1:ns1,1:ns1,1:(n.occasions-1)], dq2F[1:ns1,1:ns1,1:(n.occasions-1)])

  
  # MALES
  #~~~~~~~~~~~~
  # Calculate probability of non-encounter (dq) and reshape the array for the encounter
  # probabilities
  for (t in 1:(n.occasions-1)){
    for (s in 1:ns2){
      dp1M[s,s,t] <- po1M[s,t]
      dq1M[s,s,t] <- 1-po1M[s,t]
      dp2M[s,s,t] <- po2M[s,t]
      dq2M[s,s,t] <- 1-po2M[s,t]
    } #s
    for (s in 1:(ns2-1)){
      for (m in (s+1):ns2){
        dp1M[s,m,t] <- 0
        dq1M[s,m,t] <- 0
        dp2M[s,m,t] <- 0
        dq2M[s,m,t] <- 0
      } #s
    } #m
    for (s in 2:ns2){
      for (m in 1:(s-1)){
        dp1M[s,m,t] <- 0
        dq1M[s,m,t] <- 0
        dp2M[s,m,t] <- 0
        dq2M[s,m,t] <- 0
      } #s
    } #m
  } #t
  
  # Define the multinomial likelihood
  for (t in 1:((n.occasions-1)*ns2)){
    marrM[t,1:(n.occasions *ns2-(ns2-1)),1] ~ dmulti(pi1M[t,1:(n.occasions*ns2-(ns2-1))], released1M[t])
    marrM[t,1:(n.occasions *ns2-(ns2-1)),2] ~ dmulti(pi2M[t,1:(n.occasions*ns2-(ns2-1))], released2M[t])
  }

  # Define cell probabilities of the multistate m-array
  pi1M[1:((n.occasions-1)*ns2), 1:(n.occasions*ns2-(ns2-1))] <- marrProb(ps1M[1:ns2,1:ns2,1:(n.occasions-1)], dp1M[1:ns2,1:ns2,1:(n.occasions-1)], dq1M[1:ns2,1:ns2,1:(n.occasions-1)])
  pi2M[1:((n.occasions-1)*ns2), 1:(n.occasions*ns2-(ns2-1))] <- marrProb(ps2M[1:ns2,1:ns2,1:(n.occasions-1)], dp2M[1:ns2,1:ns2,1:(n.occasions-1)], dq2M[1:ns2,1:ns2,1:(n.occasions-1)])

})




# Bundle data
#          Year 14   15   16    17    18    19    20    21    22    23   24   # YEAR
# REPRODUCTION
R        <-  c(  0,   2,   3,    2,    6,   14,   15,   18,   19,   31,  27)  # reproductive females
J        <-  c(  0,   6,  10,    7,   22,   42,   51,   46,   61,   96,  70)  # Number of kittens
RsF      <-  c(  3,   3,   5,    4,    2,    5,    4,    3,    5,    5,   4)  # females released
RsM      <-  c(  4,   2,   4,    4,    4,    3,    6,    4,    5,    5,   4)  # males released
# COUNTS
# FEMALES
C.FRep   <-  c(  0,   2,   3,    2,    6,   14,   15,   18,   19,   31,  27)
# MALES
C.Mrep   <-  c(  0,   3,   3,    7,   13,   18,   23,   30,   43,   46,  47)  # Potentially reproductive males
# SEX-RATIO
# Kittens    
Jm       <-  c(  0,   4,   3,    3,    9,   20,   20,   20,   29,   34)       
Jf       <-  c(  0,   2,   3,    4,   11,   17,   26,   20,   24,   31) 

# Tf     <-  c(  3,   5,  10,   14,   18,   33,   42,   56,   63,   75,  86)  # Total females
# Tm     <-  c(  4,	  4,   9,	15,	  20,	30,	  47,	49,	  65,	73,	 69)  # Total males

CimM     <-  c(  0,   0,   1,    1,    1,    0,    0,    0,    0,    0,   0)
CimF     <-  c(  0,   0,   0,    0,    0,    1,    0,    0,    0,    0,   0)

mean.C <- 45
K <- 15

# DATA
str(data    <-     list(marrF=marrF,
                        marrM=marrM,
                        C.FRep = C.FRep,
						C.Mrep = C.Mrep ,
                        mean.C = mean.C,					
                        J=J,
						Jm=sum(Jm),
						Jf=sum(Jf),
                        CimM=CimM,
                        CimF=CimF,
                        R1=57,
						R2=16,
                        RsF=c(RsF, rep(5,K)),
                        RsM=c(RsM, rep(5,K))))


# CONSTANTS
# Number of true states
nts1 <- 12 # Females
nts2 <- 10 # Males
str(constants <-   list(n.occasions=ncol(CH1), 
                        ns1=nts1,
                        ns2=nts2,						
                        R=R,
                        K=15,					
                        released1F=rowSums(marrF[,,1]), 
                        released2F=rowSums(marrF[,,2]),
                        released1M=rowSums(marrM[,,1]), 
                        released2M=rowSums(marrM[,,2])))


# INITIAL VALUES
set.seed(1961)
str(Inits   <-     list(beta          =   runif(1,-2.0,2.0),
                        mean.mr1F     =   runif(1,0.01,1.0),
                        mean.mr2F     =   runif(1,0.01,1.0),
                        mean.mr3F     =   runif(1,0.01,1.0),
                        mean.mr1M     =   runif(1,0.01,1.0),                        
                        mean.mr2M     =   runif(1,0.01,1.0),
                        mean.mr3M     =   runif(1,0.01,1.0),
                        mean.mo1F     =   runif(1,0.01,1.0),
                        mean.mo2F     =   runif(1,0.01,1.0),
                        mean.mo3F     =   runif(1,0.01,1.0),
                        mean.mo1M     =   runif(1,0.01,1.0),
                        mean.mo2M     =   runif(1,0.01,1.0),
                        mean.mo3M     =   runif(1,0.01,1.0),
                        mean.FF       = c(runif(1,0.80,1.0),
                                          runif(1,0.60,1.0)),
                        mean.FM       = c(runif(1,0.80,1.0),
                                          runif(1,0.60,1.0)),										  
                        mean.p        = c(runif(1,0.90,1.0),
                                          runif(1,0.90,1.0)),
                        mean.rF       = c(runif(1,0.20,0.5),
                                          runif(1,0.20,0.9)),
                        mean.rM       = c(runif(1,0.20,0.5),
                                          runif(1,0.20,0.9)),
                        mean.alpha    =   runif(1,0.00,0.6),
                        mean.nonPhiM  = c(runif(1,0.00,0.7),
                                          runif(1,0.00,0.7)),
                        mean.nonPhiF  = c(runif(1,0.00,0.7),
                                          runif(1,0.00,0.7)),
                        mean.delta    =   runif(1,0.70,0.9),
                        mean.fec      =   runif(1,0.00,3.5),
                        log.IF        =   rnorm(26),
                        log.IM        =   rnorm(26),						
                        sigma.mr1F    =   runif(1,0.00,2.9),
                        sigma.mr2F    =   runif(1,0.00,2.3),
                        sigma.mr3F    =   runif(1,0.00,4.0),                       
                        sigma.mr1M    =   runif(1,0.00,2.3),						
                        sigma.mr2M    =   runif(1,0.00,2.3),
                        sigma.mr3M    =   runif(1,0.00,4.0),
                        sigma.mo1F    =   runif(1,0.00,3.0),
                        sigma.mo2F    =   runif(1,0.00,2.0),
                        sigma.mo3F    =   runif(1,0.00,3.0),                        
                        sigma.mo1M    =   runif(1,0.00,2.0),						
                        sigma.mo2M    =   runif(1,0.00,2.0),
                        sigma.mo3M    =   runif(1,0.00,3.0),
                        sigma.delta   =   runif(1,0.00,0.5),
                        sigma.fec     =   runif(1,0.00,0.3),
                        sigma.p       = c(runif(1,0.00,2.3),
                                          runif(1,0.00,1.7)),
                        sigma.y1      =   runif(1,0.50,3.0),
                        sigma.y2      =   runif(1,0.50,3.0)))

# PARAMETERS MONITORED
params <- c('s1F','rk1F','o1F','mean.mr1F','mean.mo1F','sigma.mr1F','sigma.mo1F',
            's1M','rk1M','o1M','mean.mr1M','mean.mo1M','sigma.mr1M','sigma.mo1M',
            's2F','s3F','FF','alpha','nonPhiF',
            's2M','s3M','FM','alpha','nonPhiM',
            'rk2F','rk3F','o2F','o3F','beta',
            'rk2M','rk3M','o2M','o3M',
            'mean.mr2F','mean.mr3F','mean.mo2F','mean.mo3F',
            'mean.mr2M','mean.mr3M','mean.mo2M','mean.mo3M',
            'delta','fec','mean.alpha','mean.fec','xi','di',
            'mean.FF','mean.FM',
            'mean.rF','mean.rM',
            'mean.delta','mean.nonPhiF','mean.nonPhiM',
            'NtH','NtM','NT','FRep',
            'N1','N2','F2','N2T','B3','N1M','N2M','N3M','IF','IM','extinct',
            'sigma.mr2F','sigma.mr3F','sigma.mo2F','sigma.mo3F',
            'sigma.mr2M','sigma.mr3M','sigma.mo2M','sigma.mo3M',
            'sigma.delta','sigma.fec',
            'sigma.y1','sigma.y2','N.nonPhiwb','N.nonPhiRs','N.Phiwb','N.PhiRs',
            'f2im','f2wb','f2rs','n2im','n2wb','n2rs','n1wb')

Rmodel <- nimbleModel(code=Code, constants=constants, data=data, inits=Inits, check=FALSE, calculate=FALSE)

# Function to initialize complex nodes
InitNod<-function(simNodes){
  simNodeScalar <- Rmodel$expandNodeNames(simNodes)
  allNodes <- Rmodel$getNodeNames()
  nodesSorted <- allNodes[allNodes %in% simNodeScalar]
  set.seed(1960) # to fix simulations
  for(n in nodesSorted) {
    Rmodel$simulate(n)
    depNodes <- Rmodel$getDependencies(n)
    Rmodel$calculate(depNodes)
  }
}

InitNod(simNodes = 'IF')
Rmodel$IF
InitNod(simNodes = 'IM')
Rmodel$IM

Cmodel <- compileNimble(Rmodel)
conf<-configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE, thin=5)
		  
MCMC <- buildMCMC(conf)
Cmcmc <- compileNimble(MCMC, project = Rmodel)

# MCMC settings
nb = 50000
ni = 50000 + nb
nc = 3
start.time2<-Sys.time()
outNim <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, #inits=Inits,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE,
                  WAIC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # running time
