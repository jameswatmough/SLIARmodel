
# two-group SLAIR model

# source the ode solvers
library("deSolve")
# ggplot2 provides routines for fancy plotting
library("ggplot2")
# tidyverse provides routines for manipulating data frames into something ggplot likes
library("tidyverse")


# ODE model

age_names = c( 'younger','middle','older' )
num_ages = length(age_names)
age_dist = c(.3,.4,.3)

status_names = c( 'unvaccinated','vaccinated' )
num_status = length(status_names)
vac_levels = c(.5,.75,.8)

population = 800000*diag(age_dist)%*%cbind(1-vac_levels,vac_levels)
dimnames(population) = list(age_names,status_names)

progression_base = 2/7
symptom_base = progression_base*(2/sqrt(1-.3) - 1)
hospital_base = progression_base*(2/sqrt(1-.2) - 1) 
Ro = 5
infectivity_base = progression_base*Ro/(1+.5*progression_base/(symptom_base+progression_base))

## State variables
# base parameters, with descriptions
# most parameters are structured by age and status, 
# some could be structured by stage as well
# structure is named in most cases for readability
# note no checking is done to ensure parameters have the correct number of values!
param = list(
	num_latent = 2,
	prog_latent = 1,        # 2 days on average with latent infection, Erlang distribution
	num_infectious = (num_infectious = 2),
	prog_infectious = progression_base,  # 7 days on average infectious, Erlang distribution
	dev_symptoms = # rate infections developing "symptoms which can't be ignored", based on 70% symptomatic 
		symptom_base*array( 
			c( 1/6,2/3,1,       # increasing with age
				1/12,1/3,1/2),     # dereasing with vaccination status
			dim=c(num_ages,num_status)
			),   
	prog_hospital = # hospitalation rates given isolation-worthy symptoms (base 20% cases hospitalized)
		hospital_base*array(
			c(1e-3, 1e-2, 1,   # status 1, stage 1
	      1e-4, 1e-4, 1,   # status 2, stage 1
	      1e-2, 1e-1, 2,   # status 1, stage 2
	      1e-3, 1e-3, 1),  # status 2, stage 2
			dim=c(num_ages,num_status,num_infectious),
			dimnames=list(age_names,status_names,NULL)
			),
	survival = 1,    # don't kill anyone
	inf_stage = infectivity_base*c(1,.5),     # second infectious stage has lower viral load, ergo less infectious
  
	inf_group = array(
		1,        # vaccine and age effect of infectiousness given infection, assumed minimal
		dim=c(num_ages,num_status),
		dimnames=list(age_names,status_names)
		),      
	sus_group =  array(
		c( 1, 1, 1,    # status 1  (unvaccinated, no reduction)
			.2,.2,.2),   # status 2  (vaccinated 80% effective)
		dim=c(num_ages,num_status),
		dimnames=list(age_names,status_names)
		),      
	contact_rate = c(1,1,1),           # equal contact rates for each group
	contact_pref = c(.6,.6,.6),          # fraction contacts within group 
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

	with(as.list(param),{
		# mixing populations assumed to be S, L, A, and R;  I and H assumed in relative isolation from main population
		# assume two groups and 2 stages in each of L, I, and A
		# caution!  the two group 2 stage assumption is hard-coded in this routine
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
		return(force_of_infection)
	})
}

force.of.infection.age.vac <- function(x,param) {

	with(as.list(param),{
		# mixing populations assumed to be S, L, A, and R; also assumed to be defined globally
		N = rowSums(x[,,mixing_stages],dims=2)

		# total number of contacts distributed across age groups 
		 C = sum(contact_rate*(1-contact_pref)*rowSums(N))


		# total number of contacts made within each group
		# Ca = contact_rate*contact_pref*rowSums(N)
		# sum relative infectivity by infectious stage for each group and then sum over status (rows)
		# note the `apply` computes dot products along the 3rd dimension (infection stage)
		# and the rowSums sums along status
		infectivity = rowSums(inf_group*apply(x[,,infectious_stages],1:2,function(x) x%*%inf_stage))
		force_of_infection = 
			contact_rate*(
				contact_pref*infectivity/rowSums(N) +                                       # within group mixing
			  (1-contact_pref)*sum(contact_rate*(1-contact_pref)*infectivity)/C  # between group mixing
				)*sus_group
		# should return a dim(num_ages,num_status) array of forces of infection

		return(force_of_infection)
	})
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

    dx[,,'S'] = -force.of.infection.age.vac(x,param)*x[,,'S']
		dx[,,'L1'] = -dx[,,'S'] + importation - prog_latent*x[,,'L1']   # latent stage 1
		dx[,,'L2'] = prog_latent*x[,,'L1'] - prog_latent*x[,,'L2']   # latent stage 2
		dx[,,'A1'] = prog_latent*x[,,'L2'] - ( prog_infectious + dev_symptoms )*x[,,'A1'] # infectious stage 1
		dx[,,'A2'] = prog_infectious*x[,,'A1'] - ( prog_infectious + dev_symptoms )*x[,,'A2']   # infectious stage 2
		dx[,,'Ra'] = prog_infectious*x[,,'A2'] # recovered from mild infection
		dx[,,'I1'] = dev_symptoms*x[,,'A1'] - (prog_infectious+prog_hospital[,,1])*x[,,'I1'] # isolated stage 1
		dx[,,'I2'] = dev_symptoms*x[,,'A2'] + prog_infectious*x[,,'I1'] - (prog_infectious+prog_hospital[,,2])*x[,,'I2'] # isolated stage 2
		dx[,,'Ri'] = prog_infectious*x[,,'I2']  # recovered from isolation 
		dx[,,'H1'] = prog_hospital[,,1]*x[,,'I1'] - prog_infectious*x[,,'H1'] # hospital stage 1
		dx[,,'H2'] = prog_hospital[,,2]*x[,,'I2'] + prog_infectious*x[,,'H1'] - prog_infectious*x[,,'H2'] # hospital stage 2
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
	# return(ode_reshape_wide(x))
	return(x)
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
			res = rbind(res,x[,c(1, 1 + i + (j-1)*num_ages + num_ages*num_status*(0:(num_stages-1)) )])
		}
	}
	res = as.data.frame(res)
	names(res) = c('time',stage_names)
	res$age=rep(age_names,each=dim(x)[1],times=num_status)
	res$status=rep(status_names,each=num_ages*dim(x)[1])
	return(res)
}	

