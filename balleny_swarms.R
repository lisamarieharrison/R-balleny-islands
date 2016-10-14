#compares krill swarms within and outside of the balleny islands
#date: 14/10/2016
#author: Lisa-Marie Harrison


if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "~/Lisa/phd/Balleny Islands/remote data/krill aggregations")
} else {
  setwd(dir = "~/phd/southern ocean/Balleny Islands/remote data/krill aggregations")
}

library(chron)

island_swarms <- read.csv("AllKrillSwarms.csv", header = T)
post_swarms   <- read.csv("AllKrillSwarmsPostBalleny.csv", header = T)

par(mfrow = c(1, 2))
hist(island_swarms$volDengPerm3)
hist(post_swarms$volDengPerm3)
hist(island_swarms$Corrected_length)
hist(post_swarms$Corrected_length)

#K-S test for swarm density
ks.test(island_swarms$volDengPerm3, post_swarms$volDengPerm3)

#create datetime
island_swarms$datetime <- chron(dates. = paste0(substr(island_swarms$Date_S, 1, 4), "/", substr(island_swarms$Date_S, 5, 6), "/", substr(island_swarms$Date_S, 7, 8)), times. = island_swarms$Time_S, format = c(dates = "y/m/d", times = "h:m:s"), out.format = c(dates = "d/m/y", times = "h:m:s"))
post_swarms$datetime <- chron(dates. = paste0(substr(post_swarms$Date_S, 1, 4), "/", substr(post_swarms$Date_S, 5, 6), "/", substr(post_swarms$Date_S, 7, 8)), times. = post_swarms$Time_S, format = c(dates = "y/m/d", times = "h:m:s"), out.format = c(dates = "d/m/y", times = "h:m:s"))

#number of hours in each area
island_hours <- as.numeric(max(island_swarms$datetime) - min(island_swarms$datetime))*24
post_hours <- as.numeric(max(post_swarms$datetime) - min(post_swarms$datetime))*24


#swarms/hour
#twice as many krill swarms per hour at Balleny Islands (may need to run by km to account for speed?)
nrow(island_swarms)/island_hours
nrow(post_swarms)/post_hours






