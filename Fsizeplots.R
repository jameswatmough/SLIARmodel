

# Run GillespieSSA for the general SLIAR model
# generate final size plots

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
  Ro = 4.0
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

# run through 6 beta values and 30 simulations each
FSize = NULL
for (b in 1:6) {
#  param$beta = rep((b+3)/18*beta,2*param$n)
	p = param
  p$phi = c(1.3*(1-b/6),.5,1.3*b/6)
# 30 simulations for 360 days
	A = repsim(
		    SIR23$S0,
				rates=SIR23$rates,
				nu=SIR23$nu,
				param=p,
				tf=360,
				SIR23$simName,
				runs=30
	)
	f = fsize(A)
	f$b = b
	f$Ro = RoSLmAI3R(p) 
	FSize = rbind(FSize,f)
}

Ro = rep(0,6)
for (b in 1:6) {
#  param$beta = rep((b+3)/18*beta,2*param$n)
	p = param
  p$phi = c(1.3*(1-b/6),.5,1.3*b/6)
  Ro[b] = RoSLmAI3R(p)
}

# plot final sizes
dev.new()
plot(
	jitter(1+FSize$RA+FSize$RI,factor=.3),
	jitter(FSize$b),
	log="x",yaxt='n',
  main = "Final Size with changing phi, N = 10,000",
  ylab = "Ro",
  xlab = "Final Size, A+I")
axis(2,at=1:6,label=as.character(round(Ro,2)))

# now set beta back to default and increase tracing
param$beta = c(rep(beta,param$n),rep(0,param$n))
FSizephi2 = NULL
for (b in 1:6) {
  param$phi = c(.5*(1-b/6),.25,.5*b/6)

# 30 simulations for 360 days
	A = repsim(
		    SIR23$S0,
				rates=SIR23$rates,
				nu=SIR23$nu,
				param=param,
				tf=360,
				SIR23$simName,
				runs=30
	)
	f = fsize(A)
	f$b = b
	FSizephi2 = rbind(FSizephi2,f)

}
