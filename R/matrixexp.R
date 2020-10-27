# compute $u\exp(Vt)v$
# where $u$ is a vector of subclinical compartments
# $V$ is the transmission matrix
# $v$ is a vector of imporation compartments
# $t$ is a quarantine time
# the product $\exp(Vt)v$ denotes the release-state of individuals initial in states $v$ at start of quarantine
# the product $u\exp(Vt)v$ denotes the fraction of those individuals still in subclinical infectious stages
# the matrix below assumes an S->L->L->(IAR)_2 progression 
library(expm)  
# or use Matrix::expm below instead,
# but note expm library is faster and more accurate

g=2/5.71
p=.2
a=2/10
n=2/10
pL1=.25
pL2=.25
pA1=.25
pA2=.25


V <- matrix(c(
	-g,      0, 0, 0, 0, 0, 0, 0,
	 g,     -g, 0, 0, 0, 0, 0, 0,
	 0,(1-p)*g,-a, 0, 0, 0, 0, 0,
	 0,      0, a,-a, 0, 0, 0, 0,
	 0,    p*g, 0, 0,-n, 0, 0, 0,
	 0,      0, 0, 0, n,-n, 0, 0,
	 0,      0, 0, a, 0, 0, 0, 0,
	 0,      0, 0, 0, 0, n, 0, 0),
 nrow=8,ncol=8,byrow=TRUE)

importation=matrix(c(pL1,pL2,0,0,pA1,pA2,0,0),ncol = 1,byrow=FALSE)
subclinical=matrix(c(1,1,0,0,1,1,0,0),nrow = 1,byrow=TRUE);

escape_prob = function(t) {subclinical%*%expm(V*t)%*%importation}
t = seq(0,15,by=.1)
plot(t,sapply(t,escape_prob),ylim=c(0,1))


