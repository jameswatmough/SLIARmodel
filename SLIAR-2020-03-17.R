# example 1

# Run GillespieSSA with the SLIAR model and compare with the solution from the ODE
# 

library("ggplot2")

source("SLIAR-SSA.R")


# random parameter values 
parms <- c(beta=.01, epsilon=0, delta=0.2, kappa=1/4, p=0.8, alpha=1/7, eta = 1/4, f = .8)

# Initial state vector for ssa
# A = 1 for the initia 'invasion'
# I = 20, A = 40 for the 'switch on distancing' scenario
S0  <- c(S=1000, L=0, I = 0, A = 1, R = 0)               
#S0  <- c(S=1000, L=0, I = 20, A = 40, R = 0)               

# final time for SSA
tf <- 60  
# output times for ode solver
times <- seq(0, tf, by=0.01)



#set.seed(1)   # use set.seed to repeat the simulation with the same random number sequence
# Run the SSA simulator 
# set verbose=FALSE to omit timing output
res.s1 <- ssa(x0 = S0,a=rates,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)

# change the transmission to reflect distancing
parms['beta'] = .004
# run 20 simulations and store the results in a big matrix
# set verbose FALSE to suppress all the output from ssa
A = NULL
for (run in seq(1,20)) {
  cat(c('starting run ',run,'\n'))
  res_ssa <- ssa(x0 = S0,
a=rates,nu,parms,tf,
method=ssa.otl(),simName,
verbose=FALSE,consoleInterval=1)
# output of ssa includes $data, which has t, S, W for each event
# append a fourth column with run number to ssa$data
# bind result to A (add it as additional rows to bottom)
  A = rbind(A,cbind(res_ssa$data,'run'=as.factor(run)))
}
# The result is a big matrix with all the outputs 
# we'll convert this to a data.frame, because data-frames are cool.
A = data.frame(A)

# use ggplot2 to make a nice plot
# probably need to change the colour scheme
ggplot(A,aes(x=t,y=I))+geom_line(aes(colour=run))

# run the ode simulator 
# one trick to get numerical solutions of population models without problems
# from round-off errors introducint negative populations is to transform the model to log
# populations,
# This means the initial conditions for the ode need to be log transformed
# initial state vector for ode (log)
x0 <- log(S0)
names(x0) = c("S","L","I","A","R")
res.o1 <- ode(y = x0, times, func = SLIARrhs, parms)

# need to figure out how to use ggplot to add the ode results in
# for now just make a new plot
dev.new()
plot(res_ode[,1],res_ode[,'I'],col='red')

