# load the libraries
# the gillespie algorithms
# first attempt to update to SSA2
require("GillespieSSA2")
# the ode solvers
require("deSolve")

SLmAIn.sim = function(m,n) {

	# The SLm(AI)nR model combines 
	# a Erlang latent distribution 
	# symptomatic and asymptomatic infectious periods
	# m latent stages
	# n infectious stages

	simName <- paste("SL",m,"(AI)",n," epidemic model", sep='')

	# compartment names to use for GillespieSSA algorithm
	Lnames = paste("L",1:m,sep="")
	Anames = paste("A",1:n,sep="")
	Inames = paste("I",1:n,sep="")

	# sample initial conditions, 
	S0 <- c(10000,1,rep(0,m-1+2*n+2))
	names(S0) = c("S",Lnames,Anames,Inames,"RA","RI")

	# one possible stochastic simulation (SSA) model for this uses 
	# five events (or transitions)
	# the model consists of specifying 
	# the changes in the state variables (transitions) 
	# and the rates for each event
	# last column is imporation of latent infections
	rates <- c(
		# incidence term based on proportional incidence
		# beta should be a vector for infectivities of A and I compartments
		paste(
			"S*sum(beta*c(",
			paste(c(Anames,Inames),collapse=","),
			"))/(",
			paste(names(S0),collapse="+"),
			")",
			sep=""
			),
		# progression through L,A and I at rate alpha with m, n and n stages
		paste("alpha*L",1:m,sep=""),
		paste("alpha*A",1:n,sep=""),
		paste("alpha*I",1:n,sep=""),
		# A->I progression at rates phi_1, ... , phi_n
		paste("phi[",1:n,"]*A",1:n,sep=""),
		# importation into L1
		"gamma"
	)   


	# The State-change matrix
	# has one column for each transition (2+m+3n) and one row for each state variable (3+m+2n)
	# last column is imporation of latent infections
	nu = matrix(
		0, 
		nrow=3+m+2*n, 
		ncol=2+m+3*n,
		dimnames = list(names(S0),NULL)
	)
	# S->L->...A->RA
	diag(nu[1:(1+m+2*n),])=-1
	diag(nu[2:(1+m+2*n),])=1
	# move the ones for An and In down to recovered row
	nu[c('I1','RA'),1+m+n]=c(0,1)
	nu['RI',1+m+2*n]=1
	# tracing (phi)
	nu[(1+m)+(1:n),(1+m+2*n)+(1:n)]=-diag(n)
	nu[(1+m+n)+(1:n),(1+m+2*n)+(1:n)]=diag(n)
	# importation (gamma)
	nu['L1',2+m+3*n]=1

	return(list(rates = rates,nu = nu,S0 = S0,simName = simName))

}

# run 20 simulations and store the results in a big matrix
# set verbose FALSE to suppress all the output from ssa

repsim = function(S0,rates,nu,param, tf, simName, runs, method = ssa.d()) {
	A = data.frame(NULL)
	for (run in seq(1,runs)) {
	res_ssa <- GillespieSSA2::ssa(
	  initial_state = S0,
	  reactions = rates,
	  nu,
	  params = param,
	  final_time = tf,
	  method=method,simName,
	  verbose=FALSE)
	# ssa.otl for tau leaping, dt for 
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

RoSLmAI3R = function(param) {
	sum( with(param, 
		   beta*c(
				   1/(phi[1]+alpha),
					 alpha/(phi[1]+alpha)/(phi[2]+alpha),
					 alpha^2/(phi[1]+alpha)/(phi[2]+alpha)/(phi[3]+alpha),
				   phi[1]/(phi[1]+alpha)/alpha,
				   ( phi[1]/(phi[1]+alpha) + phi[2]*alpha/(phi[1]+alpha)/(phi[2]+alpha) )/alpha,
				   ( phi[1]/(phi[1]+alpha) + phi[2]*alpha/(phi[1]+alpha)/(phi[2]+alpha) + phi[3]*alpha^2/(phi[1]+alpha)/(phi[2]+alpha)/(phi[3]+alpha) )/alpha)
			 ) )
}


