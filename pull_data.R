# libraries
# scripts to pull data from repositories
# scripts for various analyses and plots of data

library(ggplot2)

# read in cases from John Hopkins database
data_Global = read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

# read in NY Times data
url_nytimes = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/"
data_US = read.csv(paste(url_nytimes,"us.csv",sep=''))
data_US_state = read.csv(paste(url_nytimes,"us-states.csv",sep=''))
data_US_county = read.csv(paste(url_nytimes,"us-counties.csv",sep=''))

# read in cases from Covid19Canada database
data_Canada = read.csv("https://raw.githubusercontent.com/ccodwg/Covid19Canada/master/timeseries_prov/cases_timeseries_prov.csv")

# John Hopkins time series data are in wide format with dates as X%m.%d.%y for column names 
# the US and global files have different column numbers,  need to fix the next line to generalize away from the '5'
days = grep("^X",names(data_Global))


# recast the reporting dates for Covid19Canada as Dates
data_Canada$date_report = as.Date(data_Canada$date_report,format='%d-%m-%Y')

prov_Atlantic = c("Nova Scotia",
		  "Newfoundland and Labrador",
		  "New Brunswick",
		  "Prince Edward Island")

# pick out the columns for Canada and remove cruise ship and recovered
CAN_col = grep('Canada',data_Global$Country.Region)
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

Europe_col = data_Global$Country.Region%in%Europe_names & data_Global$Province.State==''

plot_Country = function(country,x=data_Global,d=days) {
  i = grep(country,x$Country.Region)
  plot(date,x[i,-(1:4)])
  plot(date,x[i,d],log='y')
}

plot_Prov = function(prov_name,x=data_Global,d=days) {
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
long_data = function(rows, col = 'Province.State',x = data_Global, d = days) {
  long = data.frame(NULL)
  for (prov in rows) {
		df = data.frame(
			date = as.Date(names(x)[d],format='X%m.%d.%y'), 
		  cum.cases = as.numeric(x[prov,d]),
		  prov = as.character(x[prov,col])
		)
		df$incidence = c(0,diff(df$cum.cases))
		df$mavg = filter(df$incidence,rep(1/7,7))
		long = rbind(long,subset(df,cum.cases>0))
  }
  return(long)
}

# plot given rows from the John Hopkins data sets
# use col for names for legend
plot_rows = function(rows, col = 'Province.State', x = data_Global, d = days) {
  long = long_data(rows,col = col, x=x,d=d)
  p = ggplot(long,aes(x=date,y=cum.cases))+geom_step(aes(color=prov)) + ylab("Cumulative Cases")
  print(p + coord_trans(y='log10') + scale_color_manual(values=rainbow(length(rows))))
}

plot_US_states = function(states, data = data_US_state) {
	p = ggplot(subset(data,state%in%states),aes(x=date,y=cases,col=state))
	p = p + geom_point() +coord_trans(y='log10') 
	p = p + geom_smooth(method='lm',formula=y~x)
	return(p)
}

plot_Province_daily = function(prov = c("New Brunswick")) {
	ggplot(subset(data_Canada,province%in%prov),aes(x=date_report,y=cases)) + geom_step(aes(color=province)) + ylab("Daily Cases")
}

doubling.time = function(t,y,window) { 
	# compute doubling times from vectors of date (t) and incidence (y)
	n = length(window)
	# result will be a vector of length(y)-n
	dtime = rep(0,length=length(y)-n)
  for (i in 1:length(dtime)) {
		Y = log(y[i+window],2)
		X = as.integer(t[i+window])
		dtime[i] = (n*sum(X*Y)-sum(Y)*sum(X))/(n*sum(X^2) - sum(X)^2) 
	}
	return(data.frame(date = t[-(window)],dtime = dtime))
}
