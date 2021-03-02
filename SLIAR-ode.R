# A minimalist script to demonstrate deSolve with an SIR-type model


# source the ode solvers
library("deSolve")


library("ggplot2")

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

paramset <- read.csv("data/parameter-sets.csv")

param = paramset[paramset$scenario=='Diamond',]
# output times for ode solver
tf <- 60  
times <- seq(0, tf, by=0.01)
S0  <- c(S=1000, L=0, I = 1, A = 0, R = 0)               

# run the ode simulator 
res_ode <- ode(y = S0, times, func = SLIAR.ode, param)

plot(res_ode[,1],res_ode[,'I'],col='red')
