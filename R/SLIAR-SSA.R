
# Second, set up the ode solver for the same model with the same parameters

# load the libraries
# the gillespie algorithms
library("GillespieSSA")
# the ode solvers
library("deSolve")

# this has been superseded by SLmAIn.init 
# and can be replaced with 
#   sim = SLmAIn.init(m=1,n=1)
# the plots in the libraries make use of a simulation name
SLIAR.sim = function(N0=10000) {

	simName <- "SLIAR epidemic model"

	# sample initial conditions, 
	S0 <- c(N0-1,1,rep(0,m-1+2*n+2))

	# The SLIAR model combines 
	# a latent period 
	# symptomatic and asymptomatic infectious periods
	# a removed compartment

	# one possible stochastic simulation (SSA) model for this uses 
	# five events (or transitions)
	# the model consists of specifying 
	# the changes in the state variables (transitions) 
	# and the rates for each event
	# last column is imporation of latent infections
	rates <- c(
		"beta*S*(epsilon*L + I + delta*A)/(S+L+I+A+R)",    
		"p*kappa*L",
		"(1-p)*kappa*L",
		"eta*A",
		"f*alpha*I",
		"(1-f)*alpha*I",
		"gamma"
	)   


	# The State-change matrix
	# has one column for each transition and one row for each state variable
	# last column is imporation of latent infections
	nu  <- matrix(
		c(-1,  0, 0,   0,  0,  0, 0,
      +1, -1, -1,  0,  0,  0,+1,
       0, +1,  0,  0, -1, -1, 0,
       0,  0, +1, -1,  0,  0, 0,
       0,  0,  0, +1, +1,  0, 0),
    nrow=5,byrow=TRUE
	) 
	return(list(rates = rates,nu = nu,S0 = S0,simName = simName))
}

# ODE model
# set up the right hand side for the ode solver
# using `with` is the prefered method, but allowing 
# for missing parameters allows use of global variables
# note the switch to log variables to avoid numerical problems near zero populations
SLIAR.ode <- function(t,x, parms=NULL) {
	with(as.list(c(x,parms)),{

		dS = -beta*S*(epsilon*L + I + delta*A)/(S+L+I+A+R);
		dL = gamma + beta*S*(epsilon*L + I + delta*A)/(S+L+I+A+R) - kappa*L;
		dI = p*kappa*L - alpha*I;
		dA = (1-p)*kappa*L - eta*A;
		dR = f*alpha*I + eta*A;
		list(c(dS,dL,dI,dA,dR))

	})

}


SLIAR.Ro = function(param) {
	Ro = with(param, beta*(epsilon/kappa + p/alpha + (1-p)*delta/eta))
  names(Ro) = param$scenario
	return(Ro)
}


# run 20 simulations and store the results in a big matrix
# set verbose FALSE to suppress all the output from ssa

repsim = function(S0,rates,nu,param, tf, simName, runs, method = ssa.d()) {
	A = data.frame(NULL)
	for (run in seq(1,runs)) {
	res_ssa <- ssa(x0 = S0,
		a=rates,nu,parms = param,tf,
		method=method,simName,
		verbose=FALSE)
	# ssa.otl for taul leaping, dt for 
# output of ssa includes $data, which has t, S, L, ... for each event
# append a fourth column with run number to ssa$data
# bind result to A (add it as additional rows to bottom)
		A = rbind(A,data.frame(res_ssa$data,run=as.factor(run)))
	}
	return(A)
}

pluckrun = function(data, sim, col = 3:5) {
	wide = subset(data, run==sim)
	long = stack(wide[,col])
	long$t = rep(wide$t,times=length(col))
	return(long)
}

fsize = function(A) {
	# find the index of the final time for each sim
	# this assumes the times are in increasing order
	tails = rep(0,length(unique(A$run)))
	for (i in as.integer(unique(A$run))) { 
		tails[i] = sum(A$run==i)
	}
	# the indices are relative to each subset
	# use a cumulative sum to find the indices for A
	return(A[Reduce('+',tails,accumulate=TRUE),])
}


