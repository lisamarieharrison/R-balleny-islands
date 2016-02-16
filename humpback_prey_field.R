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
library(Matching) #K-S bootstrap
library(plotrix) #vectorField

#source required functions
function_list <- c("gcdHF.R",
                   "deg2rad.R"
)

for (f in function_list) {
  source(paste("C:/Users/43439535/Documents/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/", f, sep = ""))
}

#subset sightings to only HB whales (HB = Species 07)
sighting <- subset(sighting, Species == 07)

#subset sightings to only MI platform
sighting <- subset(sighting, Platform == "MI")

#plot of sightings by date
plot(table(sighting$Date), ylab = "Number of sightings", main = "HB sightings by date") 

sighting$Date <- chron(dates. = as.character(sighting$Date), format = "d/m/y")
sighting$Time <- chron(times. = as.character(sighting$Time), format = "h:m:s")

effort$Date <- chron(dates. = as.character(effort$Date), format = "d/m/y")
effort$Time <- chron(times. = as.character(effort$Time), format = "h:m:s")

#subset sightings and effort to only krill dates
sighting <- sighting[chron(dates. = "3/2/2015", format = "d/m/y") <= sighting$Date & sighting$Date <= chron(dates. = "6/2/2015", format = "d/m/y"), ]
effort   <- effort[chron(dates. = "3/2/2015", format = "d/m/y") <= effort$Date & effort$Date <= chron(dates. = "6/2/2015", format = "d/m/y"), ]

#remove error values
krill$arealDen[krill$arealDen == -9.900000e+37] <- NA

#effort where MI observers present
#next zero cell included to give total time on effort
#doesn't include changes in observers 
#doesn't include half effort for now
#effort <- effort[effort$NObserversM != c(tail(effort$NObserversM, -1), 2), ]
effort$NObserversM[effort$NObserversM == 1] <- 0 #to count half effort as off effort
effort$NObserversM[effort$EffortStatus == "OF"] <- 0 #remove observers if waypoint is marked as being off transect
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

#remove duplicated krill rows
krill <- krill[!duplicated(krill), ]

#format krill times and dates to match sightings
krill$Ping_date <- chron(dates. = as.character(krill$Ping_date), format = "y-m-d", out.format = "d/m/y")
krill$Ping_time <- chron(times. = as.character(krill$Ping_time), format = "h:m:s")

krill$datetime    <- chron(dates. = krill$Ping_date, times. = krill$Ping_time, format = c(dates = "d/m/y", times = "h:m:s"))
sighting$datetime <- chron(dates. = sighting$Date, times. = sighting$Time, format = c(dates = "d/m/y", times = "h:m:s"))


#remove krill when off effort for whales
krill_on_effort <- NULL
krill_datetime_on_effort <- NULL
krill_lat_on_effort <- NULL
krill_long_on_effort <- NULL
for (i in 1:length(start_datetime)) {
  
  w <- start_datetime[i] < krill$datetime & end_datetime[i] > krill$datetime
  
  if (sum(w) > 0) {
    krill_on_effort <- c(krill_on_effort, krill$arealDen[w])
    krill_datetime_on_effort <- c(krill_datetime_on_effort, as.character(krill$datetime[w]))
    krill_lat_on_effort <- c(krill_lat_on_effort, krill$Latitude[w])
    krill_long_on_effort <- c(krill_long_on_effort, krill$Longitude[w])
  }
  
}
krill_datetime_on_effort <- chron(dates. = substr(krill_datetime_on_effort, start = 2, stop = 9), times. = substr(krill_datetime_on_effort, start = 11, stop = 18), format = c(dates = "d/m/y", times = "h:m:s"))

plot(krill_datetime_on_effort, krill_on_effort, pch = 19, xlab = "Date", ylab = "krill density gm2")
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
whale_present <- rep(FALSE, length(krill_datetime_on_effort))
whale_number  <- rep(0, length(krill_datetime_on_effort))
sighting_id   <- NULL
for (i in 1:nrow(sighting)) {
  w <- which.min(abs(sighting$datetime[i] - krill_datetime_on_effort)*24*60)
  if (abs(sighting$datetime[i] - krill_datetime_on_effort[w])*24*60 <= 15) {
    whale_present[w] <- TRUE
    whale_number[w]  <- sighting$BestNumber[i]
    sighting_id[w]   <- sighting$Index[i]
  }
}

#plot krill density at cells with and without whales
boxplot(krill_on_effort ~ whale_present, ylab = "krill density gm2")
legend("topright", c("TRUE = whale sighting in cell", "FALSE = no whale in cell"), bty = "n")

boxplot(log(krill_on_effort) ~ whale_present, ylab = "log(krill density gm2)")
legend("topright", c("TRUE = whale sighting in cell", "FALSE = no whale in cell"), bty = "n")

#two sample t-test
#do cells with more whales have a higher krill density?
#log transformed because of skew
krill_on_effort[krill_on_effort == 0] <- NA #remove single 0 value
t.test(log(krill_on_effort[whale_present]), log(krill_on_effort[!whale_present]), alternative = "greater")

#Kolmogorov-Smirnov Test for non-parametric data
#do cells with whales present have not less than (at least equal to) the krill density where whales are absent?
#note that specifying alternative is done in opposite way to t.test
ks.test(krill_on_effort[whale_present], krill_on_effort[!whale_present], alternative = "less")
#use K-S bootstrap implementation because there are ties present
ks.boot(krill_on_effort[whale_present], krill_on_effort[!whale_present], nboots = 1000, alternative = "less")

#------------------------ PLOT KRILL DENSITY AGAINST NUMBER OF WHALES ----------------------#

#number of whales varies between 1 and 6
plot(krill_datetime_on_effort, krill_on_effort, pch = 19, xlab = "Date", ylab = "krill density gm2", xaxt = "n")
axis(1, sighting$datetime, sighting$BestNumber, col.ticks = "red") #ticks at whale locations with number
title("Krill density with number of whales in each sighting in red")
legend("topright", col = "red", "Whale sighting location", lwd = 2, bty= "n")

plot(krill_on_effort, whale_number, pch = 19, xlab = "Krill density gm2", ylab = "Number of whales in sighting", cex.lab = 2)
points(na.omit(krill_on_effort[sighting_id == na.omit(sighting$Index[sighting$Behaviour == 5])]), na.omit(whale_number[sighting_id == na.omit(sighting$Index[sighting$Behaviour == 5])]), col = "red")
legend("topright", col = c("red", "black"), pch = 19, c("Feeding", "Other behaviour"), bty = "n")


#------------------------------------ PREDICTIVE GLM --------------------------------------#

krill.glm <- glm(whale_present ~ krill_on_effort, family = binomial(link = logit))
summary(krill.glm)

table(whale_present[!is.na(krill_on_effort)], round(krill.glm$fitted.values))

#------------------------------------ SHIP HEADING ---------------------------------------#

#gps doesn't take into account on effort times

direction <- gps$Heading[gps$PCTime %in% c("3/02/2015", "4/02/2015", "5/02/2015", "6/02/2015")]
x <- gps$Longitude[gps$PCTime %in% c("3/02/2015", "4/02/2015", "5/02/2015", "6/02/2015")]
y <- gps$Latitude[gps$PCTime %in% c("3/02/2015", "4/02/2015", "5/02/2015", "6/02/2015")]

plot(krill_long_on_effort, krill_lat_on_effort)
vectorField(direction, 1, x, y, scale = 0.005, vecspec = "deg")
points(krill_long_on_effort, krill_lat_on_effort, col = "red", pch = 19)



