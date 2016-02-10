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
