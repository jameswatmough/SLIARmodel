


# Run GillespieSSA for the general SLIAR model
# for comparison to Gaia's 30-day 'September Ends' scenarios

library("ggplot2")
library("GillespieSSA")

source("R/SLmAIn-SSA.R")

# set up very simple parameter set 
# 2 latent and 3 infectious stages
param = list(
  n = 3, 
  m = 2,
  # mean time from infection to recovery
  T = 12,
# set some tracing and no imporation
  phi = c(.4,.3,.2),
  gamma = 0,
  Ro = 3.1
)

# choose alpha based on infection period and numbers of stages
param$alpha = with(param, (n+m)/T)
# with this setup, Ro, with phi=0, is n*beta/alpha = beta*T*n/(n+m)
# so set beta based on Ro

# scale beta to specified Ro assuming phi=0 (no tracing, all cases undetected)
beta = with(param, Ro*(n+m)/n/T)

param$beta = c(rep(beta,param$n),rep(0,param$n))
#param$beta = rep(.5*beta,2*param$n)


SIR23 = SLmAIninit(param$m,param$n)

# insert the model for quarantine
# multipy beta by 1/2 for I1+I2+I3>6
# and by 1/3 for I1+I2+I3>20
SIR23$rates[1] = "ifelse(I1+I2+I3>5,ifelse(I1+I2+I3<20,0.5,1/3),1)*S*sum(beta*c(A1,A2,A3,I1,I2,I3))/(S+L1+L2+A1+A2+A3+I1+I2+I3+RA+RI)"

p = param

SIR23$S0[1] = 20000
LastObs = NULL
runs = 200
for (delay in 4:6) {
	SIR23$S0[2:6] = 0
	SIR23$S0[delay] = 1

	A = repsim(
		    SIR23$S0,
				rates=SIR23$rates,
				nu=SIR23$nu,
				param=p,
				tf=30,
				SIR23$simName,
				runs=runs
	)

	for (k in 1:runs) { 
		B = subset(A,run==k)
		LastObs = rbind(LastObs,cbind(B[dim(B)[1],],delay=delay))
	}
}
LastObs$SubClinical = rowSums(LastObs[,5:7])
LastObs$Clinical = rowSums(LastObs[,8:10])
LastObs$Latent = rowSums(LastObs[,3:4])
LastObs$Recovered = rowSums(LastObs[,11:12])
LastObs$Total = rowSums(LastObs[,c(8:10,12)])


p = ggplot(LastObs,aes(y = Total, x=c('Stage 1','Stage 2','Stage 3')[7-delay])) 
p+ geom_boxplot(outlier.shape=NA) + geom_jitter(width=.2,height=.2) + labs(y="Cumulative Cases at 30 days",x="") 
