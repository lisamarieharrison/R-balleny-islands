#balleny islands krill prey field analysis
#author: Lisa-Marie Harrison
#date: 10/02/2016
#3 passes over the same transect so time and date used to match sightings rather than lat and long

setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/csv")
gps      <- read.csv("GpsData.csv", header = T)
sighting <- read.csv("Sighting.csv", header = T)
env      <- read.csv("Environment.csv", header = T)
effort   <- read.csv("Effort.csv", header = T)
krill    <- read.csv("Krill.csv", header = T)
library(chron)
library(ggplot2)
library(Matching) #ks.boot
library(plotrix) #vectorField
library(geosphere) #destPoint

#source required functions
function_list <- c("gcdHF.R",
                   "deg2rad.R"
)

for (f in function_list) {
  source(paste("C:/Users/43439535/Documents/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/", f, sep = ""))
}

#subset sightings to only HB whales (HB = Species 07) & to only MI platform
sighting <- subset(sighting, Species == 07 & Platform == "MI")

#plot of sightings by date
plot(table(sighting$Date), ylab = "Number of sightings", main = "HB sightings by date") 

#remove error values
krill$arealDen[krill$arealDen == -9.900000e+37] <- NA

#remove duplicated krill rows
krill <- krill[!duplicated(krill), ]

#format krill times and dates to match sightings
krill$datetime    <- chron(dates. = as.character(krill$Ping_date), times. = as.character(krill$Ping_time), format = c(dates = "y-m-d", times = "h:m:s"), out.format = c(dates = "d/m/y", times = "h:m:s"))
sighting$datetime <- chron(dates. = as.character(sighting$Date), times. = as.character(sighting$Time), format = c(dates = "d/m/y", times = "h:m:s"))
effort$datetime   <- chron(dates. = as.character(effort$Date), times. = as.character(effort$Time), format = c(dates = "d/m/y", times = "h:m:s"))
gps$datetime      <- chron(dates. = as.character(gps$PCTime), times. = gps$Time, format = c(dates = "d/m/y", times = "h:m:s"))

sighting    <- subset(sighting, sighting$datetime >= min(krill$datetime) & sighting$datetime <= max(krill$datetime))
gps         <- subset(gps, gps$datetime >= min(krill$datetime) & gps$datetime <= max(krill$datetime))
effort      <- subset(effort, effort$datetime >= min(krill$datetime) & effort$datetime <= max(krill$datetime))

#effort where MI observers present
#next zero cell included to give total time on effort
#doesn't include changes in observers or half effort
effort$NObserversM[effort$NObserversM == 1 | effort$EffortStatus == "OF"] <- 0 #to count half effort as off effort
w <- which(effort$NObserversM == 0 & c(0, head(effort$NObserversM, -1)) == 2 | effort$NObserversM == 2 & c(0, head(effort$NObserversM, -1)) == 0)
effort <- effort[w, ]
effort <- effort[-nrow(effort), ]

#gps index for start and end of effort
gps_full <- effort$GpsIndex[effort$NObserversM == 2] #start of full effort
gps_end  <- effort$GpsIndex[effort$NObserversM == 0] #end of effort

#time and date of start and end of effort

start_datetime <- chron(dates. = as.character(gps$PCTime[na.omit(match(gps_full, gps$Index))]), times. = as.character(gps$Time[na.omit(match(gps_full, gps$Index))]), format = c(dates = "d/m/y", times = "h:m:s"))
end_datetime   <- chron(dates. = as.character(gps$PCTime[na.omit(match(gps_end, gps$Index))]), times. = as.character(gps$Time[na.omit(match(gps_end, gps$Index))]), format = c(dates = "d/m/y", times = "h:m:s"))

#reminder to check start and end time length for duplicates
if (length(start_datetime) != length(end_datetime)) {
  stop("Effort start and end times do not match up")
}

#add 15 mins buffer to each to get full end bins
start_datetime <- start_datetime - 15/60/24 
end_datetime   <- end_datetime + 15/60/24

#functions to subset to only on effort times
withinTimes <- function(x, start, end) {
  
  #check if a row is between any on effort times
  #dataRow: row of a data frame
  #start: vector of start date and times as a chron object
  #end: vector of end date and times as a chron object
  
  datetime_chron <- x[which(names(x) == "datetime")]
  
  datetime <- chron(dates. = substr(datetime_chron, 2, 9), times. = substr(datetime_chron, 11, 18), 
                  format = c(dates = "d/m/y", times = "h:m:s"))

  return (any(start < datetime & end > datetime)) 
  
}

onEffort <- function(data, start, end) {
  
  #for every row in a data frame, check whether it is within any on effort times
  #if FALSE, remove the row, else move on to next row
  #data: data frame with each row an observation
  #start: vector of start date and times as a chron object
  #end: vector of end date and times as a chron object
  
  data_oneffort <- data[apply(data, 1, FUN = withinTimes, start = start_datetime, end = end_datetime), ]
  
  return (data_oneffort)
  
}

#subset krill and gps to only on effort times
krill <- onEffort(krill, start_datetime, end_datetime)
gps   <- onEffort(gps, start_datetime, end_datetime)


plot(krill$datetime, krill$arealDen, pch = 19, xlab = "Date", ylab = "krill density gm2")
rug(sighting$datetime, ticksize = 0.03, side = 1, lwd = 0.5, col = "red", quiet = TRUE) #ticks at whale locations
title("Krill density with whale sightings in red")
legend("topright", col = "red", "Whale sighting location", lwd = 2, bty = "n")

#plot sightings along transect
plot(krill$Longitude, krill$Latitude, type = "l", xlab = "Longitude", ylab = "Latitude")
points(gps$Longitude[gps$Index %in% sighting$GpsIndex], gps$Latitude[gps$Index %in% sighting$GpsIndex], col = "red", pch = 19)

d <- qplot(Longitude, Latitude, data=krill, colour= arealDen)
d + scale_colour_gradient(low = "grey", high = "blue", name = "Krill density gm2") + 
  theme_bw() +
  geom_point(data = gps[gps$Index %in% sighting$GpsIndex & gps$Latitude < -66, ], aes(Longitude, Latitude), colour = "red", shape = 8)



#------------------------ TEST KRILL DENSITY WITH AND WITHOUT WHALES ----------------------#

#does a krill cell contain a whale?
#for each sighting, which cell is the whale closest to?
whale_present <- rep(FALSE, nrow(krill))
whale_number  <- rep(0, nrow(krill))
sighting_id   <- NULL
for (i in 1:nrow(sighting)) {
  w <- which.min(abs(sighting$datetime[i] - krill$datetime)*24*60)
  if (abs(sighting$datetime[i] - krill$datetime[w])*24*60 <= 15) {
    whale_present[w] <- TRUE
    whale_number[w]  <- sighting$BestNumber[i]
    sighting_id[w]   <- sighting$Index[i]
  }
}

#plot krill density at cells with and without whales
boxplot(krill$arealDen ~ whale_present, ylab = "krill density gm2")
legend("topright", c("TRUE = whale sighting in cell", "FALSE = no whale in cell"), bty = "n")

boxplot(log(krill$arealDen) ~ whale_present, ylab = "log(krill density gm2)")
legend("topright", c("TRUE = whale sighting in cell", "FALSE = no whale in cell"), bty = "n")

#two sample t-test
#do cells with more whales have a higher krill density?
#log transformed because of skew
krill$arealDen[krill$arealDen == 0] <- NA #remove single 0 value
t.test(log(krill$arealDen[whale_present]), log(krill$arealDen[!whale_present]), alternative = "greater")

#Kolmogorov-Smirnov Test for non-parametric data
#do cells with whales present have not less than (at least equal to) the krill density where whales are absent?
#note that specifying alternative is done in opposite way to t.test
ks.test(krill$arealDen[whale_present], krill$arealDen[!whale_present], alternative = "less")
#use K-S bootstrap implementation because there are ties present
ks.boot(krill$arealDen[whale_present], krill$arealDen[!whale_present], nboots = 1000, alternative = "less")

#------------------------ PLOT KRILL DENSITY AGAINST NUMBER OF WHALES ----------------------#

#number of whales varies between 1 and 6
plot(krill$datetime, krill$arealDen, pch = 19, xlab = "Date", ylab = "krill density gm2", xaxt = "n")
axis(1, sighting$datetime, sighting$BestNumber, col.ticks = "red") #ticks at whale locations with number
title("Krill density with number of whales in each sighting in red")
legend("topright", col = "red", "Whale sighting location", lwd = 2, bty= "n")

plot(krill$arealDen, whale_number, pch = 19, xlab = "Krill density gm2", ylab = "Number of whales in sighting", cex.lab = 2)
points(na.omit(krill$arealDen[sighting_id == na.omit(sighting$Index[sighting$Behaviour == 5])]), na.omit(whale_number[sighting_id == na.omit(sighting$Index[sighting$Behaviour == 5])]), col = "red")
legend("topright", col = c("red", "black"), pch = 19, c("Feeding", "Other behaviour"), bty = "n")


#------------------------------------ PREDICTIVE GLM --------------------------------------#

krill.glm <- glm(whale_present ~ krill$arealDen, family = binomial(link = logit))
summary(krill.glm)

table(whale_present[!is.na(krill$arealDen)], round(krill.glm$fitted.values))

#------------------------------------ SHIP HEADING ---------------------------------------#

direction <- gps$Heading
x <- gps$Longitude
y <- gps$Latitude

plot(krill$Longitude, krill$Latitude)
vectorField(direction, 1, x, y, scale = 0.005, vecspec = "deg")
points(krill$Longitude, krill$Latitude, col = "red", pch = 19)


#find true lat and long of sighting
#estimate monkey island height of 25m and 0.03 degrees per 0.1 reticle
sighting$Reticles[sighting$Reticles == 0] <- NA
reticle_distance <- 25*tan(deg2rad(90 - sighting$Reticles/0.1*0.03))


#if all sightings were 5km away where would they be on the map?


sightingAngle <- function(x, gps) {
  
  #calculate the absolute heading of a sighting taking into account ships heading
  #x = row of sighting
  #gps = full gps matrix
  
  angle <- as.numeric(x[which(names(x) == "Angle")])
  index <- as.numeric(x[which(names(x) == "GpsIndex")])
  
  if (length(gps$Heading[gps$Index == index]) > 0) {
    
    #lhs
    if (angle <= 90) {
      angle_true <- angle + gps$Heading[gps$Index == index]
    } else {
      angle_true <- gps$Heading[gps$Index == index] - (360 - angle)
    }
    
    if(angle_true < 0) {
      angle_true <- 360 + angle_true
    }
    
  } else {
    return (NA)
  }
  
  return (angle_true)
  
}


sighting$angle_true <- apply(sighting, 1, sightingAngle, gps = gps)


sightingLatLong <- function (x, gps, distance) {
  
  #calculates the true latitude and longitude of an object from its bearing, distance and point of observation
  #x = row of sighting
  #gps = full gps matrix
  #distance = distance to object in m
  
  index <- as.numeric(x[which(names(x) == "GpsIndex")])
  angle <- as.numeric(x[which(names(x) == "angle_true")])
  
  if (length(gps$Longitude[gps$Index == index]) > 0) {
    
    lon <- destPoint(p = c(gps$Longitude[gps$Index == index], gps$Latitude[gps$Index == index]),
                     b = angle, d = distance)[1]
    
    lat <- destPoint(p = c(gps$Longitude[gps$Index == index], gps$Latitude[gps$Index == index]),
                     b = angle, d = distance)[2]
  } else {
    
    lon <- NA
    lat <- NA
    
  }
  
  return(c(lat, lon))
  
}

true_lat_long <- data.frame(t(apply(sighting, 1, sightingLatLong, gps = gps, distance = 5000)))
colnames(true_lat_long) <- c("lat", "long")


plot(krill$Longitude, krill$Latitude, pch = 19, xlab = "Longitude", ylab = "Latitude")
points(gps$Longitude[gps$Index %in% sighting$GpsIndex], gps$Latitude[gps$Index %in% sighting$GpsIndex], col = "red", pch = 19)
points(true_lat_long$lon, true_lat_long$lat, col = "orange", pch = 19)








