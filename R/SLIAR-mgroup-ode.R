
# two-group SLAIR model

# source the ode solvers
library("deSolve")
# ggplot2 provides routines for fancy plotting
library("ggplot2")
# tidyverse provides routines for manipulating data frames into something ggplot likes
library("tidyverse")


# ODE model

age_names = c( 'younger','older' )
num_ages = length(age_names)

status_names = c( 'unvaccinated','vaccinated' )
num_status = length(status_names)

population = 800000*matrix(
	c(.4*(1-.6),.4*.6,
		.6*(1-.8),.6*.8),
	nrow = 2,
	byrow=TRUE,
	dimnames=list(age_names,status_names)
	)

## State variables
# base parameters, with descriptions
param = list(
	num_latent = 2,
	prog_latent = 1,        # 2 days on average with latent infection, Erlang distribution
	num_infectious = 2,
	prog_infectious = 2/7,  # 7 days on average infectious, Erlang distribution
	dev_symptoms = 2/7*c(1/3,1),   # 25% and 50% of cases develop "symptoms which can't be ignored", need to separate out vaccinated
	prog_hospital_1 = 2/7*c(1/19,1/4),  # hospitalation rates (based on fractions hospitalized)
	prog_hospital_2 = 2/7*c(1/19,1/4),  # hospitalation rates (based on fractions hospitalized)
	survival = 1,
	inf_stage = c(1,.5),     # second infectious stage has lower viral load, ergo less infectious
  
	inf_group = array(
		1,        # vaccine effect of infectiousness given infection 
		dim=c(num_ages,num_status),
		dimnames=list(age_names,status_names)
		),      
	sus_group =  array(
		c(.55,.28,.55,.28),      # random numbers for now!
		dim=c(num_ages,num_status),
		dimnames=list(age_names,status_names)
		),      
	contact_rate = c(1,1),           # equal contact rates for each group
	contact_pref = c(.6,.6),          # fraction contacts within group 
	importation = array(
		0,        # vaccine effect of infectiousness given infection 
		dim=c(num_ages,num_status),
		dimnames=list(age_names,status_names)
		)      
)	

stage_names = c(
 'S',
 'L1', # latent group-stage
 'L2',
 'A1', # mild (a)symptomatic and infectious group-stage
 'A2',
 'Ra', # recovered with mild symptoms
 'I1', # symptomatic and presumed isolating group-stage
 'I2', 
 'Ri', # recovered from isolation 
 'H1', # hospitalized, group-stage
 'H2',
 'Rh'  # recovered from hospital 
)

num_stages = length(stage_names)
mixing_stages = c(1:6,9,12)
infectious_stages = 4:5

ics = array(
	data = 0.,
  dim = c(num_ages,num_status,num_stages),
  dimnames = list(age_names,status_names,stage_names)
	)
ics[,,'S'] = population
ics['younger','unvaccinated','L1'] = 1


force.of.infection.age <- function(x,param) {

		# mixing populations assumed to be S, L, A, and R;  I and H assumed in relative isolation from main population
		# assume two groups and 2 stages in each of L, I, and A
		# caution!  the two group 2 stage assuption is hard-coded in this routine
		N = c(
			sum(x[mixing_stages]),
		  sum(x[num_stages+mixing_stages])
			)
		# total number of contacts made across both groups 
		C = sum(contact_rate*(1-contact_pref)*N)
		infectivity = c(
			sum(inf_stage*x[infectious_stages]),
			sum(inf_stage*x[num_stages+infectious_stages])
			)
		force_of_infection = 
			contact_rate*contact_pref*sus_group*inf_group*infectivity/N +                                       # within group mixing
			contact_rate*(1-contact_pref)*sus_group*sum(contact_rate*(1-contact_pref)*inf_group*infectivity)/C  # between group mixing
}

force.of.infection.age.vac <- function(x,param) {

		# mixing populations assumed to be S, L, A, and R; also assumed to be defined globally
		N = rowSums(x[,,mixing_stages],dims=2)


		# total number of contacts distributed across both groups 
		C = sum(contact_rate*(1-contact_pref)*rowSums(N))


		# total number of contacts made within each group
		Ca = contact_rate*contact_pref*rowSums(N)
		# sum relative infectivity by infectious stage for each group
		infectivity = rowSums(inf_group*rowSums(inf_stage*x[,,infectious_stages],dims=2))
		force_of_infection = 
			contact_rate*(
				contact_pref*infectivity/rowSums(N) +                                       # within group mixing
			  (1-contact_pref)*sum(contact_rate*(1-contact_pref)*infectivity)/C  # between group mixing
				)*sus_group
		# should return a dim(num_ages,num_status) array of forces of infection

		return(force_of_infection)
}

# set up the right hand side for the ode solver
# using `with` is the prefered method, but allowing 
# for missing parameters allows use of global variables
mgroup.ode <- function(
	t,                     # current time
	x,                     # named vector of current states
	parms=NULL)            # named vector of parameters
{
	# apply indexing to the state vector for imporved readability
	dim(x) = c(num_ages,num_status,num_stages)
	dimnames(x) = list(age_names,status_names,stage_names) 

	dx = array(
		0,
  	dim = c(num_ages,num_status,num_stages), 
  	dimnames = list(age_names,status_names,stage_names) 
  )

	# unwind states and parameters to make code a bit more readable
	# and return right hand side of ode as a named list
	with(as.list(parms),{

    dx[,,'S'] = -force.of.infection.age.vac(x,param)*x[c(1,13)]
		dx[,,'L1'] = -dx[,,'S'] + importation - prog_latent*x[,,'L1']   # latent stage 1
		dx[,,'L2'] = prog_latent*x[,,'L1'] - prog_latent*x[,,'L2']   # latent stage 2
		dx[,,'A1'] = prog_latent*x[,,'L2'] - ( prog_infectious + dev_symptoms )*x[,,'A1'] # infectious stage 1
		dx[,,'A2'] = prog_infectious*x[,,'A1'] - ( prog_infectious + dev_symptoms )*x[,,'A2']   # infectious stage 2
		dx[,,'Ra'] = prog_infectious*x[,,'A2'] # recovered from mild infection
		dx[,,'I1'] = dev_symptoms*x[,,'A1'] - (prog_infectious+prog_hospital_1)*x[,,'I1'] # isolated stage 1
		dx[,,'I2'] = dev_symptoms*x[,,'A2'] + prog_infectious*x[,,'I1'] - (prog_infectious+prog_hospital_2)*x[,,'I2'] # isolated stage 2
		dx[,,'Ri'] = prog_infectious*x[,,'I2']  # recovered from isolation 
		dx[,,'H1'] = prog_hospital_1*x[,,'I1'] - prog_infectious*x[,,'H1'] # hospital stage 1
		dx[,,'H2'] = prog_hospital_2*x[,,'I2'] + prog_infectious*x[,,'H1'] - prog_infectious*x[,,'H2'] # hospital stage 2
		dx[,,'Rh'] = survival*prog_infectious*x[,,'H2'] # recovered from hospital 

	return(list(as.vector(dx)))

	}) # close with 
	# flatten back to state vector and return
} # close function


# deSolve ode routines need initial states and a list of output times

tf <- 360  
times <- seq(0, tf, length=120)

# run the ode simulator and return the results in a wide format grouped by age and status
ode_sim = function(
	param,
	y = ics,    # vector of initial states
	times = seq(0,tf,length=120),         # vector of output times
	func = mgroup.ode
) {
	# run default ode solver
	x = ode(y, times, func, param) 
	# rewrap the results in long format
	return(ode_reshape_wide(x))
}

ode_reshape_long = function(res) {
	# really just a reminder of how to reshape the wide version of the results
	pivot_longer(res,all_of(stage_names),names_to='stage')
}

ode_reshape_wide = function(x) {
	# reshape to data frame grouped by age and status and columns for each infection stage
	res = NULL
	for (i in 1:num_ages) {
		for (j in 1:num_status) {
			res = rbind(res,x[,c(1,1+j+num_ages*num_status*(0:(num_stages-1)))])
		}
	}
	res = as.data.frame(res)
	names(res) = c('time',stage_names)
	res$age=rep(age_names,each=dim(x)[1],times=num_status)
	res$status=rep(status_names,each=num_ages*dim(x)[1])
	return(res)
}	

