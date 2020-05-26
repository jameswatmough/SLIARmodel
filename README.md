# Simple SIR-type ode and branching process simulations

#  Files

 * *SLIAR-testruns.R* sources a model and parameter set, runs the ode and SSA routines, and produces a few simple plots

 * *SLIAR-importation.R* sources a model and parameter set, sets Ro below 1 and sets importation rate to 3 latent infections per day, runs the SSA routines, and produces a few simple plots

 * *R/SLIAR-SSA.R* model definitions and base scripts for SSA and ode runs

 * *data/parameter-sets.csv* sample parameter sets

# data structures

  * *param* row from parameter-sets for current computations
	* *rates* list of transition rates for SSA routine
	* *nu*  matrix of transitions, one column for each transition, one row for each state variable
	* *rhs* rates for the ode approximation
