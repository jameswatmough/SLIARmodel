# create a scatter plot of deaths by date and zone indicating number of deaths

library(ggplot2)

deaths = read.csv('data/deaths-nb.csv')
deaths$date = as.Date(deaths$date)

ggplot(deaths,aes(x=date,y=as.factor(zone))) + geom_count() + labs(x="Date",y="Zone",title="CoViD Deaths in NB")

