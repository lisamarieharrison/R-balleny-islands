#balleny islands krill prey field analysis
#author: Lisa-Marie Harrison
#date: 10/02/2016

setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/csv")

gps      <- read.csv("GpsData.csv", header = T)
sighting <- read.csv("Sighting.csv", header = T)
env      <- read.csv("Environment.csv", header = T)
effort   <- read.csv("Effort.csv", header = T)

#subset sightings to only HB whales (HB = Species 07)
sighting <- subset(sighting, Species == 07)

#plot of sightings by date
plot(table(sighting$Time), ylab = "Number of sightings", main = "HB sightings by date")

#effort where MI observers present
#next zero cell included to give total time on effort
#doesn't include changes in observers
effort <- effort[effort$NObserversM != c(tail(effort$NObserversM, -1), 2), ]
effort <- effort[-1, ]

gps_full <- effort$GpsIndex[effort$NObserversM == 2] #start of full effort
gps_half <- effort$GpsIndex[effort$NObserversM == 1] #start of half effort
gps_end  <- effort$GpsIndex[effort$NObserversM == 0] #end of effort
