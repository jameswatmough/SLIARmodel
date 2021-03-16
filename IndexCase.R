# sample script to generate states at time of first observed case

# Run GillespieSSA for the general SLIAR model
# generate distributions of infections by stage at time of first observed case

library("ggplot2")
library("GillespieSSA")
library("dplyr")

source("R/SLmAIn-SSA.R")

# set up very simple parameter set 
# 2 latent and 3 infectious stages
param = list(
	S0 = 10000,
  n = 3, 
  m = 2,
  # mean time from infection to recovery
  T = 12,
# set some tracing and no imporation
  phi = c(.2,.4,.6),
  gamma = 0,
  Ro = 4.0,
	runs = 2000
)

# choose alpha based on infection period and numbers of stages
param$alpha = with(param, (n+m)/T)
# with this setup, Ro, with phi=0, is n*beta/alpha = beta*T*n/(n+m)
# so set beta based on Ro

# scale beta to specified Ro assuming phi=0 (no tracing, all cases undetected)
beta = with(param, Ro*(n+m)/n/T)

param$beta = c(rep(beta,param$n),rep(0,param$n))
#param$beta = rep(.5*beta,2*param$n)


SIR23 = SLmAIn.init(m = param$m, n=param$n)

A = repsim(
		    SIR23$S0,
				rates=SIR23$rates,
				nu=SIR23$nu,
				param=param,
				tf=60,
				SIR23$simName,
				runs=param$runs
)

# subset the row from each run with the first observed case
FirstObs = NULL
for (k in 1:param$runs) { 
  B = subset(A,run==k)
	# look for first row with something in the I columns (8:10)
  FirstObs = rbind(FirstObs,B[match(1,rowSums(B[,8:10])),])
}
# collapse the substages into single counts
FirstObs$Infectious = as.integer(rowSums(FirstObs[,5:7]))
FirstObs$Latent = as.integer(rowSums(FirstObs[,3:4]))
FirstObs$Recovered = as.integer(rowSums(FirstObs[,11:12]))
FirstObs$Infected = as.integer(rowSums(FirstObs[,3:7]))

count_first_obs = count(FirstObs,Infectious,Latent,name='density')
count_first_obs$density = count_first_obs$density/param$runs
# 
p = ggplot(count_first_obs,aes(x=Infectious,y=Latent,fill=density)) 
p = p + geom_tile() 
p = p + scale_fill_gradient(trans="log",breaks = 5*10^(-(4:1)),labels= 5*10^(-(4:1)))
p = p + geom_text(aes(label=density))

p + scale_y_continuous(
    breaks = function(x) unique(
      floor(pretty(seq(0, (max(x) + 1) * 1.1)))))

# and plot the total number of infected (L+I)
dev.new()
 ggplot(count(FirstObs,Infected),aes(x=Infected,y=n))+geom_step()

