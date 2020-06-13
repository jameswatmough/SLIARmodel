# libraries

library(ggplot2)

# read in cases from John Hopkins database
datapath = "~/projects/SARS-CoV-2/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
setwd(datapath)
system("git pull")
filename = paste(datapath,"time_series_covid19_confirmed_global.csv",sep="")
data = read.csv(filename)

datapathCA = "~/projects/SARS-CoV-2/data/Covid19Canada/timeseries_prov/"
filenameCAprov =paste(datapathCA,"cases_timeseries_prov.csv",sep="")
dataCA = read.csv(filenameCAprov)

dataCA$date_report = as.Date(dataCA$date_report,format='%d-%m-%Y')
days = 5:dim(data)[2]

date = as.Date(names(data)[-(1:4)],format='X%m.%d.%y')

prov_Atlantic = c("Nova Scotia",
		  "Newfoundland and Labrador",
		  "New Brunswick",
		  "Prince Edward Island")

# pick out the columns for Canada and remove cruise ship and recovered
CAN_col = grep('Canada',data$Country.Region)
CAN_col = CAN_col[-c(3,12:13)]
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

Europe_col = data$Country.Region%in%Europe_names & data$Province.State==''

plot_Country = function(country,x=data,d=days) {
  i = grep(country,x$Country.Region)
  plot(date,x[i,-(1:4)])
  plot(date,x[i,d],log='y')
}

plot_Prov = function(prov_name,x=data,d=days) {
  i = grep(prov_name,x$Province.State)
  prov_data = data.frame(date = as.Date(names(x)[d],format='X%m.%d.%y'), 
		  cum.cases = as.numeric(x[i,d]))
  p = ggplot(subset(prov_data,cum.cases>0),aes(x=date,y=cum.cases))+
      geom_point()+
      xlab("Date") + ylab("Cumulative Cases") +
      coord_trans(y='log10') 
  return(p)
}

# pull out canadian rows and put in long format
# generalize this to strip out a few rows and convert them to long format
long_data = function(rows, col = 'Province.State',x = data, d = days) {
  long = data.frame(NULL)
  for (prov in rows) {
  df = data.frame(date = as.Date(names(x)[d],format='X%m.%d.%y'), 
		  cum.cases = as.numeric(x[prov,d]),
		  prov = as.character(x[prov,col]))
  long = rbind(long,subset(df,cum.cases>0))
  }
  return(long)
}

plot_rows = function(rows, x = data, d = days) {
  long = long_data(rows,x=x,d=d)
  p = ggplot(long,aes(x=date,y=cum.cases))+geom_step(aes(color=prov)) + ylab("Cumulative Cases")
  print(p + coord_trans(y='log10') + scale_color_manual(values=rainbow(length(rows))))
}

plot_NB = function() {
	ggplot(subset(dataCA,province=="New Brunswick"),aes(x=date_report,y=cases)) + geom_step(aes(color=province))
}
