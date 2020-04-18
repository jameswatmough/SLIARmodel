# libraries

library(ggplot2)
library(viridis)

# read in cases from John Hopkins database
datapath = "~/projects/SARS-CoV-2/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
filename = paste(datapath,"time_series_covid19_confirmed_global.csv",sep="")
x = read.csv(filename)

days = dim(x)[2]

date = as.Date(names(x)[-(1:4)],format='X%m.%d.%y')

# pick out the columns for Canada and remove cruise ship and recovered
CAN_col = grep('Canada',x$Country.Region)
CAN_col = CAN_col[-c(3,12:13)]
CAN_col = grep('Canada',x$Country.Region)
Europe_names = c('Austria',
		'Albania',
		'Andorra',
		'Belarus',
		'Belgium',
		'Bulgaria',
		'Bosnia and Herzegovina',
		'Bulgaria',
		'Croatia',
		'Czechia',
		'Denmark',
		'Finland',
		'France',
		'Germany',
		'Greece',
		'Hungary',
		'Iceland',
		'Ireland',
		'Italy',
		'Latvia',
		'Liechtenstein',
		'Lithuania',
		'Luxembourg',
		'Malta',
		'Monaco',
		'Montenegro',
		'Netherlands',
		'North Macedonia',
		'Norway',
		'Poland',
		'Portugal',
		'Romania',
		'Russia',
		'San Marino',
		'Serbia',
		'Slovakia',
		'Slovenia',
		'Spain',
		'Sweden',
		'Switzerland',
		'Ukraine',
		'United Kingdom')

Europe_rows = x$Country.Region%in%Europe_names & x$Province.State==''

JPN_col = grep('Japan',x$Country.Region)

plot(date,x[140,-(1:4)])
dev.new()
plot(date,x[140,5:days],log='y')


# New Brunswick Data
nb = data.frame(date = as.Date(names(x)[5:days],format='X%m.%d.%y'), cum.cases = as.numeric(x[CAN_col[5],5:days]))
p = ggplot(subset(nb,cum.cases>0),aes(x=date,y=cum.cases))+geom_point()+xlab("Date")+ylab("Cumulative Cases") + coord_trans(y='log10') 
dev.new()
print(p)

# pull out canadian rows and put in long format
can_long = data.frame(NULL)
for (prov in CAN_col) {
  df = data.frame(date = as.Date(names(x)[5:days],format='X%m.%d.%y'), cum.cases = as.numeric(x[prov,5:days]),prov = as.character(x$Province.State[prov]))
	can_long = rbind(can_long,subset(df,cum.cases>0))
}

p = ggplot(can_long,aes(x=date,y=cum.cases))+geom_step(aes(color=prov)) + ylab("Cumulative Cases")
dev.new()
print(p + coord_trans(y='log10') + scale_color_viridis(discrete=TRUE))

# or just plot NB data
p = ggplot(subset(can_long,prov=="New Brunswick"),aes(x=date,y=cum.cases))+geom_step(aes(color=prov)) + ylab("Cumulative Cases")
p = p + coord_trans(y='log10')
print(p)

prov_Atlantic = c("Nova Scotia","Newfoundland and Labrador","New Brunswick","Prince Edward Island")
p = ggplot(subset(can_long,prov%in%prov_Atlantic),aes(x=date,y=cum.cases))+geom_step(aes(color=prov)) + ylab("Cumulative Cases") +  coord_trans(y='log10')
