# quick script to produce final size results from SSA for SLIAR model
# used with arino et al importation risk paper

library("ggplot2")

source("R/SLIAR-SSA.R")

param = list(gamma=1,beta=1.2,epsilon=0.054,delta=0.23,p=.2,kappa=1/5.71,alpha=1/10,eta = 1/10,f=0.978)
param$scenario = 'Arino et al 2020'

param$Ro = SLIAR.Ro(param)

S0  <- c(S=10000, L=0, I = 0, A = 0, R = 0)               

# run through 6 beta values and 30 simulations each
FSize = NULL
p = param
for (Ro in seq(.5,2.5,by=.25) ) {
#  param$beta = rep((b+3)/18*beta,2*param$n)
  p$beta = param$beta*Ro/param$Ro
	A = repsim(
		S0,
		rates=SLIAR.rates,
		nu=SLIAR.nu,
		param=p,
		tf=92,
		simName,
		runs=50
	)
	f = fsize(A)
	f$Ro = Ro
	FSize = rbind(FSize,f)
}
