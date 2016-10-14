#compares krill swarms within and outside of the balleny islands


setwd("C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/remote data/krill aggregations")


island_swarms <- read.csv("AllKrillSwarms.csv", header = T)

k38 <- read.csv("038kHzaggregationsPostBalleny.csv", header = T)
k120 <- read.csv("120kHzaggregationsPostBalleny.csv", header = T)
