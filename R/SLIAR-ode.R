# A minimalist script to demonstrate deSolve with an SIR-type model
# age structure added 2021-08-01

# source the ode solvers
library("deSolve")
# ggplot2 provides routines for fancy plotting
library("ggplot2")
# tidyverse provides routines for manipulating data frames into something ggplot likes
library("tidyverse")

# start with the basic SLIAR model and simple examples

# ODE model
# set up the right hand side for the ode solver
# using `with` is the prefered method, but allowing 
# for missing parameters allows use of global variables
# note the switch to log variables to avoid numerical problems near zero populations
SLIAR.ode <- function(
	t,                     # current time
	x,                     # named vector of current states
	parms=NULL)            # named vector of parameters
{
	# unwind states and parameters to make code a bit more readable
	# and return right hand side of ode as a named list
	with(as.list(c(x,parms)),{    
	list(c(
		dS = -beta*S*(epsilon*L + I + delta*A)/(S+L+I+A+R),
		dL = gamma + beta*S*(epsilon*L + I + delta*A)/(S+L+I+A+R) - kappa*L,
		dI = p*kappa*L - alpha*I,
		dA = (1-p)*kappa*L - eta*A,
		dR = f*alpha*I + eta*A
	)) # close return list
	}) # close with 
} # close function

# read in parameter values from a file and pick one to play with
paramset <- read.csv("../data/parameter-sets.csv")

param = paramset[paramset$scenario=='Diamond',]

# deSolve ode routines need initial states and a list of output times

tf <- 60  
times <- seq(0, tf, by=0.01)
S0  <- c(S=1000, L=0, I = 1, A = 0, R = 0)               

# run the ode simulator and return the results in a way ggplot likes

ode_long = function(
	y = c(S=1000, L=0, I=1, A=0, R=0),    # vector of initial states
	times = seq(0,tf,length=200),         # vector of output times
	func = SLIAR.ode,
	# parameter defaults from princess diamond outbreak
	param = c(
		gamma=0,       # importation rate
		beta=1.265,    # transmission
		epsilon=0.2,   # relative infectivity of presymptomatic stage
		delta=0.2,     # relative infectivity of asymptomatic stage
		kappa=0.33,    # duration of presymptomatic stage
		p=0.25,        # proportion syptomatic
		alpha=0.33,    # duration of symptomatic stage
		eta=0.33,      # duration of asymptomatic stage
		f=0.8          # survival fraction of symptomatic
		)
	)
{
	# run default ode solver
	x = ode(y, times, func,param) 
	# rewrap the results in long format
	as.data.frame(x) %>% pivot_longer(!time,names_to='state')
}

## plot infected states vs time on single pair of axes

state_plot = function(data) {
	ggplot(data,aes(x=time,y=value)) + 
		geom_line(aes(col=state))
}

## facet plots (array of state-vs-time plots)
facet_plot = function(data) {
	ggplot(data,aes(x=time,y=value))+
		geom_line() + 
		facet_wrap(~state,scales='free')
} 

