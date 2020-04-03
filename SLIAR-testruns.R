# Run GillespieSSA with the SLIAR model and compare with the solution from the ODE
# 

library("ggplot2")

source("R/SLIAR-SSA.R")


# read in parameter values based on various scenarios

paramset <- read.csv("data/parameter-sets.csv")


Ro = with(paramset, beta*(epsilon/kappa + p/alpha + (1-p)*delta/eta) )

# use the second parameter set  (Diamond Princess paper)
param = paramset[2,2:9]

# Initial state vector for ssa
# A = 1 for the initial 'invasion'
# I = 20, A = 40 for the 'switch on distancing' scenario
S0  <- c(S=1000, L=0, I = 1, A = 0, R = 0)               
#S0  <- c(S=1000, L=0, I = 20, A = 40, R = 0)               

# final time for SSA
tf <- 60  
# output times for ode solver
times <- seq(0, tf, by=0.01)

# 20 simulations 
A = repsim(S0,rates,nu,param=param,tf=360,simName,runs=20)

# plot R for every sim
sampleruns <- ggplot(A,aes(x=t,y=R))+geom_line(aes(color = run))
print(sampleruns)

# plot L, I, A for one sim, should probably pick one that doesn't crash right away:)
B = pluckrun(A,3)
dev.new()
singlerun <- ggplot(B,aes(x=t,y=values))+geom_line(aes(color = ind))
print(singlerun)
# run the ode simulator 
# one trick to get numerical solutions of population models without problems
# from round-off errors introducint negative populations is to transform the model to log
# populations,
# This means the initial conditions for the ode need to be log transformed
# initial state vector for ode (log)
#x0 <- log(S0)
#names(x0) = c("S","L","I","A","R")
res_ode <- ode(y = S0, times, func = SLIARrhs, param)

# need to figure out how to use ggplot to add the ode results in
# for now just make a new plot
dev.new()
plot(res_ode[,1],res_ode[,'I'],col='red')

