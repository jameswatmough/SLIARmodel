# First, set up a Gillespie simulation of the SLIAR model 
# Second, set up the ode solver for the same model with the same parameters

# load the libraries
# the gillespie algorithms
library("GillespieSSA")
# the ode solvers
library("deSolve")

# the plots in the libraries make use of a simulation name
simName <- "SLIAR epidemic model"

# The SLIAR model combines 
# a latent period 
# symptomatic and asymptomatic infectious periods
# a removed compartment

# one possible stochastic simulation (SSA) model for this uses 
# five events (or transitions)
# the model consists of specifying 
# the changes in the state variables (transitions) 
# and the rates for each event
rates <- c("beta*S*(epsilon*L + I + delta*A)",    
           "p*kappa*L",
           "(1-p)*kappa*L",
           "eta*A",
           "f*alpha*I",
           "(1-f)*alpha*I")   


# The State-change matrix
# has one column for each transition and one row for each state variable
nu  <- matrix(c(-1,  0, 0,   0,  0,  0,
                +1, -1, -1,  0,  0,  0,
                 0, +1,  0,  0, -1, -1,
                 0,  0, +1, -1,  0,  0,
                 0,  0,  0, +1, +1,  0),
                 nrow=5,byrow=TRUE) 

# ODE model
# set up the right hand side for the ode solver
# using `with` is the prefered method, but allowing 
# for missing parameters allows use of global variables
# note the switch to log variables to avoid numerical problems near zero populations
SLIARrhs <- function(t,x, parms=NULL) {
  with(as.list(c(x,parms)),{

    dS = -beta*S*(epsilon*L + I + delta*A);
    dL = beta*S*(epsilon*L + I + delta*A) - kappa*L;
    dI = p*kappa*L - alpha*I;
    dA = (1-p)*kappa*L - eta*A;
    dR = f*alpha*I + eta*A;
    list(c(dS,dL,dI,dA,dR))

  })

}


