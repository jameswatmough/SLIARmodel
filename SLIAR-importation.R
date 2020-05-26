
# Run GillespieSSA for the SLIAR model with imporation 
# 

library("ggplot2")

source("R/SLIAR-SSA.R")


# read in parameter values based on various scenarios

paramset <- read.csv("data/parameter-sets.csv")


paramset$Ro = SLIAR.Ro(paramset)

# use the second parameter set  (Diamond Princess paper)
param = paramset[paramset$scenario=='Diamond',]

# set some distancing and imporation
param$beta = .5
param$gamma = 3
param$Ro = SLIAR.Ro(param)

# Initial state vector for ssa
# A = 1 for the initial 'invasion'
# I = 20, A = 40 for the 'switch on distancing' scenario
S0  <- c(S=1000, L=0, I = 1, A = 0, R = 0)               
#S0  <- c(S=1000, L=0, I = 20, A = 40, R = 0)               


# 20 simulations for 360 days
A = repsim(S0,rates=SLIAR.rates,nu=SLIAR.nu,param=param,tf=360,simName,runs=20)

# plot R for every sim
sampleruns <- ggplot(A,aes(x=t,y=R))+geom_line(aes(color = run))
print(sampleruns)

# plot L, I, A for one sim, should probably pick one that doesn't crash right away:)
B = pluckrun(A,3)
dev.new()
singlerun <- ggplot(B,aes(x=t,y=values))+geom_line(aes(color = ind))
print(singlerun)

