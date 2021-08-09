# Simple SIR-type ode and branching process simulations

## Multi-group impimentation with age and vaccine strucuture

The file 'R/SLIAR-mgroup-ode.R' contains code for a simple SLIAR model with two structures.  The defaults assume the structures are age (younger, middle, older) and vaccine status (vaccinated, unvaccinated).  Erlang distributions of sojourn times are implimented by doubling the L,I,and A compartments.  The default parameters assume L is latent, I is infectious, but isolated, and A is infectious and non-isoloting (mild, easily-dismssed symptoms).  Thus, Ra, Ri, and Rh are recovered asymptomatic, isolating, and hospitalized, respectively.  
Note the SLIAR structure is really S->L->{AIH}->R where A,I and H are parallel routes of progression and everyone begins their infectious period without symptoms (i.e., in A).

Run the model with default parameters as follows:

~~~
source('R/SLIAR-mgroup-ode.R')
res = ode_sim(param)
res_long = ode_reshape_long(res)
~~~

'res' will contain a wide version of the output, 'res_long' contains a long version

plot results with commands like
~~~
gplot(
  subset(res,status=='vaccinated'),
	aes(x=time,y=Rh)
	) +
	geom_line(aes(col=age))
~~~
