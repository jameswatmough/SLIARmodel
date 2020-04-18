# test run to think about community spread vs travel spread

# start with the imperial college parameters and the simple SLIAR model
# add in importation of presymptomatic (L) individuals at rate gamma
# look at differences in growth with Ro just above and just below 1

library("ggplot2")
source("R/SLIAR-SSA.R")

# imperial college numbers with low importation (gamma) and about 75% reduction in contacts (small beta and therefore Ro < 1)
param = c(gamma = .3,
					beta = .2,
					epsilon = .2,
					delta = .2,
					kappa = .33,
					p = .7,
					alpha = .33,
					eta = .33,
					f = .8)

# large community of 100K people
S0  <- c(S=100000, L=0, I = 1, A = 0, R = 0)               

# pick Ro and scale beta from that

Ro = 1.05
param['beta'] = with(as.list(param), Ro/(epsilon/kappa + p/alpha + (1-p)*delta/eta) )

# 20 simulations of about 100 days
A = repsim(S0,rates,nu,param=param,tf=100,simName,runs=20)

# plot R for every sim
sampleruns <- ggplot(A,aes(x=t,y=R))+geom_line(aes(color = run))
# pretty the plot up
runA <- sampleruns +  xlab("Date")+ylab("Recovered Cases") 

# if Ro > 1 with a longer run, use a log transform
# runA <- sampleruns 
#           + xlab("Date")
#           + ylab("Recovered Cases") 
#					 + coord_trans(y='log10') 
#					 + scale_y_continuous(breaks=c(20,50,100,200,500,1000,2000,5000,10000))

print(runA)  # to display the plot

# ggsave(plot = runA, file = "test-2020-04-16-Ro1.05.png")   # save the plot

# repeat with beta picked to set Ro near .9
# admire the variability between the 20 runs with the same parameters
