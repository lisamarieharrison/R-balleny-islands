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
reticle  <- read.csv("reticle.csv", header = T)
library(chron)
library(ggplot2)
library(Matching) #ks.boot
library(plotrix) #vectorField
library(geosphere) #destPoint
library(pscl) #hurdle
library(caret) #sensitivity/specificity
library(flux) #auc
library(maptools) #gcDestination
library(raster) 
library(AER) #dispersiontest
library(randomForest)
library(spgwr)

#source required functions
function_list <- c("gcdHF.R",
                   "deg2rad.R",
                   "rocCurve.R",
                   "draw_map_scale.R"
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
env$GpsTime.1     <- substr(as.POSIXct(env$GpsTime.1, format = "%I:%M:%S %p", tz = "GMT"), 12, 19)
env$datetime      <- chron(dates. = as.character(env$Time), times. = env$GpsTime.1, format = c(dates = "d/m/y", times = "h:M:S %p"))

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

#calculate time between gps readings when on effort
gps$bin_time <- rep(NA, nrow(gps))
for (i in 2:nrow(gps)) {
  bin_time <- as.numeric(gps$datetime[i] - gps$datetime[i - 1])*24*60
  
  if (bin_time > 20) {
    bin_time <- as.numeric(gps$datetime[i + 1] - gps$datetime[i])*24*60
  }
  gps$bin_time[i] <- bin_time
}

plot(krill$datetime, krill$arealDen, pch = 19, xlab = "Date", ylab = "krill density gm2")
rug(sighting$datetime, ticksize = 0.03, side = 1, lwd = 0.5, col = "red", quiet = TRUE) #ticks at whale locations
title("Krill density with whale sightings in red")
legend("topright", col = "red", "Whale sighting location", lwd = 2, bty = "n")

#plot sightings along transect
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
    whale_number[w]  <- whale_number[w] + sighting$BestNumber[i]
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

#hurdle model for count
krill.hurdle <- hurdle(whale_number ~ krill$arealDen, dist = "poisson", zero.dist = "binomial", link = "logit")
summary(krill.hurdle)
        
#------------------------------------ TRUE SIGHTING LOCATION ---------------------------------------#

#plot ship direction
direction <- gps$Heading
x <- gps$Longitude
y <- gps$Latitude

plot(krill$Longitude, krill$Latitude)
vectorField(direction, 1, x, y, scale = 0.005, vecspec = "deg")
points(krill$Longitude, krill$Latitude, col = "red", pch = 19)


#find true lat and long of sighting using reticle
#average human eye height on monkey island height is 15.08m

sightingDistance <- function(x, reticle) {
  
  #calculates true line distance to sighting using height of monkey island and reticle lookup table
  #x: row of sighting matrix for apply function
  #reticle: reticle distance lookup table
  
  measured_reticle <- as.numeric(x[which(names(x) == "Reticles")])
  
  if (is.na(measured_reticle)) {
    
    hypotenuse <- as.numeric(x[which(names(x) == "Est.Distance")])/1000
    
  } else {
    
    hypotenuse <- reticle$km[which.min(abs(reticle$angle - measured_reticle))]
    
  }
  
  distance <- sqrt(hypotenuse^2 - 0.01508^2)

  return(distance)
  
}
 
sighting$distance <- unlist(apply(sighting, 1, sightingDistance, reticle = reticle))


#sighting location given specified distance


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

sightingLatLong <- function (x, gps) {
  
  #calculates the true latitude and longitude of an object from its bearing, distance and point of observation
  #x = row of sighting
  #gps = full gps matrix
  #distance = distance to object in m
  
  index <- as.numeric(x[which(names(x) == "GpsIndex")])
  angle <- as.numeric(x[which(names(x) == "angle_true")])
  distance <- as.numeric(x[which(names(x) == "distance")])*1000
  
  if (length(gps$Longitude[gps$Index == index]) > 0) {
    
    lon <- destPoint(p = c(gps$Longitude[gps$Index == index], gps$Latitude[gps$Index == index]),
                     b = angle, d = distance)[1]
    
    lat <- destPoint(p = c(gps$Longitude[gps$Index == index], gps$Latitude[gps$Index == index]),
                     b = angle, d = distance)[2]
  } else {
    
    lon <- NA
    lat <- NA
    
  }
  
  true_lat_long <- cbind(lat, lon)
  names(true_lat_long) <- c("Latitude", "Longitude")
  
  return(true_lat_long)
  
}

sighting$angle_true <- apply(sighting, 1, sightingAngle, gps = gps)
true_lat_long <- data.frame(t(apply(sighting, 1, sightingLatLong, gps = gps)))


p <- ggplot() + geom_point(data = krill, aes(x = Longitude, y = Latitude, size = 2, colour = log(arealDen)))
p + scaleBar(lon = 165, lat = -66.3, distanceLon = 5, distanceLat = 2, distanceLegend = 5, dist.unit = "km", orientation = FALSE) + 
  #geom_point(aes(x = gps$Longitude[gps$Index %in% sighting$GpsIndex], y = gps$Latitude[gps$Index %in% sighting$GpsIndex]), color = "red") + 
  geom_point(data = true_lat_long, aes(x = Longitude, y = Latitude, size = 2), shape = 8, color = "red") + 
  theme_bw() +
  scale_color_gradient(low="blue", high="yellow", na.value="white")


#--------------------------- WEIGHTED KRILL DENSITY AROUND EACH SIGHTING -------------------------------#

timeDifference <- function (sighting, krill) {
  
  #calculates time in hours between each sighting and each krill bin
  #sighting: full sighting matrix with datetime column
  #krill: full krill matrix with datetime column
  
  diff_hours <- matrix(NA, nrow = nrow(krill), ncol = nrow(sighting))
  
  for (i in 1:length(sighting$datetime)) {
    for (j in 1:length(krill$datetime)) {
      diff_hours[j, i] <- abs(as.numeric(sighting$datetime[i] -  krill$datetime[j])*24)
    }
  }
  
  return(diff_hours)
  
}

distFromPoint <- function (x, krill, gps, truePosition=FALSE) {
  
  #returns distance in kms of each krill bin from each sighting
  #x: row of sighting from apply function
  #krill: full matrix of krill locations
  #gps: gps lookup matrix for whale sightings
  #truePosition: is the x matrix used the true position of the sighting or the observation location (default)
  #set truePosition=TRUE if using true_lat_long with colnames "Latitude" and "Longitude"
  
  distance <- NULL
  
  if(!truePosition) {
    
    origin_long <- deg2rad(gps$Longitude[gps$Index == as.numeric(x[which(names(x) == "GpsIndex")])])
    origin_lat  <- deg2rad(gps$Latitude[gps$Index == as.numeric(x[which(names(x) == "GpsIndex")])])
    
  } else {
    
    origin_long <- deg2rad(as.numeric(x[which(names(x) == "Longitude")]))
    origin_lat  <- deg2rad(as.numeric(x[which(names(x) == "Latitude")]))
    
  }
  
  for (j in 1:length(krill$Latitude)) {
    distance[j] <- gcdHF(origin_lat, origin_long, deg2rad(krill$Latitude)[j], deg2rad(krill$Longitude)[j])    
  }
  
  distance[distance > 10000] <- NA

  return(distance)
  
}

krillWeighted <- function(distance, threshold, time) {
  
  #calculates the weighted mean of krill around a sighting within a specified distance interval
  #distance: matrix of distance to each krill bin from each sighting from distFromPoint function
  #threshold: the threshold distance in km from each sighting of krill bins to use
  #time: time_difference matrix of time (hrs) between all sightings and krill bins
  
  krill_mean <- NULL
  for (i in 1:ncol(distance)) {
    weight <- 1/distance[, i]
    weight[weight < 1/threshold] <- 0
    weight[time[, i] > 1] <- 0 #exclude krill measured >1 hr since sighting
    krill_mean[i] <- weighted.mean(x = krill$arealDen, w = weight, na.rm = TRUE) 
  }
  
  return(krill_mean)
  
}

time_difference <- timeDifference(sighting, krill)

#using observation location
distance   <- apply(sighting, 1, FUN = distFromPoint, krill = krill, gps = gps)
krill_mean <- krillWeighted(distance, threshold = 5, time = time_difference)

plot(krill_mean, sighting$BestNumber, xlab = "mean krill density (gm2)", ylab = "Count in sighting")
title("Weighted mean krill density in 5km radius around sightings")

#using estimated sighting location
sighting$angle_true <- apply(sighting, 1, sightingAngle, gps = gps)
true_lat_long <- data.frame(t(apply(sighting, 1, sightingLatLong, gps = gps)))
distance   <- apply(true_lat_long, 1, FUN = distFromPoint, krill = krill, gps = gps, truePosition = TRUE)
krill_mean <- krillWeighted(distance, threshold = 5, time = time_difference)

plot(krill_mean, sighting$BestNumber, pch = 19, xlab = "mean krill density (gm2)", ylab = "Count in sighting")
title("Weighted mean krill density in 5km radius around sightings")

#poisson glm for sighting count and density
count.glm <- glm(sighting$BestNumber ~ krill_mean, family = "poisson")
summary(count.glm)


#--------------------------- IS THERE A WHALE WITHIN 15KM? -------------------------------#

#is there a whale within 5km of a krill bin?
distance[time_difference > 1] <- NA #remove distances for sightings >1hr from krill bin
whale_present <- rep(FALSE, nrow(distance))
closest_whale <- apply(distance, 1, min, na.rm = TRUE)
closest_whale[closest_whale == Inf] <- NA
whale_present[closest_whale < 3] <- TRUE

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


#glm model for presence/absence
krill.glm <- glm(whale_present ~ log(krill$arealDen), family = "binomial")
summary(krill.glm)

whale_estimate <- as.logical(round(krill.glm$fitted.values))
table(whale_present[!is.na(krill$arealDen)], whale_estimate)

sensitivity(as.factor(whale_estimate), reference = as.factor(whale_present[!is.na(krill$arealDen)]))
specificity(as.factor(whale_estimate), reference = as.factor(whale_present[!is.na(krill$arealDen)]))

#plot ROC curve
M.ROC <-rocCurve(krill.glm, threshold=0.5, data = whale_present[!is.na(krill$arealDen)]) 

par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
plot(M.ROC[1, ], M.ROC[2, ], lwd = 2, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", cex.lab = 2, cex.axis = 2)
title("ROC curve")
lines(c(0, 1), c(0, 1), col = "red")

#calculate the area under the ROC curve (0.5 = bad, 0.8 = good, 0.9 = excellent, 1 = perfect)
auc(M.ROC[1,], M.ROC[2,])
  
  
#number of individuals within 5km
countIndividuals <- function (x, sighting, threshold) {
  
  #calculates the number of individual whales within a specified distance of each krill bin
  #x: row of distance matrix for apply function
  #sighting: full sighting matrix
  #threshold: maximum distance (km) for a sighting to be included
  
  count <- sum(sighting$BestNumber[which(x < threshold)])
  
  return (count)
  
}

whale_number <- apply(distance, 1, countIndividuals, sighting = sighting, threshold = 5)

plot(krill$arealDen, whale_number, pch = 19, xlab = "krill density gm2", ylab = "number of whales")
title("Number of whales within 5km of a krill bin")

#hurdle model for count
#significant for presence/absence but not count
krill.hurdle <- hurdle(whale_number ~ log(krill$arealDen), dist = "poisson", zero.dist = "binomial", link = "logit")
summary(krill.hurdle)


#-------------------------- IS THIS CELL THE CLOSEST TO A WHALE? -----------------------------#

closest_bin <- apply(distance, 2, which.min)
closest_bin[is.na((closest_bin == Inf))] <- NA
closest_bin <- unlist(closest_bin)

closest_measured_krill <- krill$arealDen[closest_bin]

whale_present <- rep(FALSE, nrow(distance))
whale_present[closest_bin] <- TRUE

whale_number <- rep(0, nrow(distance))
whale_number[na.omit(unique(closest_bin))] <- table(closest_bin)

ks.test(krill$arealDen[whale_present], krill$arealDen[!whale_present], alternative = "less")

plot(krill$arealDen, whale_number)


closest_whale <- apply(distance, 1, min, na.rm = TRUE)

plot(log(krill$arealDen), closest_whale, pch = 19)


group <- rep(1, length(whale_present))
arealDen <- krill$arealDen
arealDen[arealDen == 0] <- NA
x <- krill$Latitude
y <- krill$Longitude

krill.glmm <- glmmPQL(whale_present ~ log(arealDen), random =~ 1|group, family = binomial, correlation =  corExp(form = ~ x + y))
summary(krill.glmm)

whale_estimate <- as.logical(round(krill.glmm$fitted.values))
table(whale_present[!is.na(arealDen)], whale_estimate)


#-------------------- PROBABILITY OF THERE BEING A SIGHTING GIVEN KRILL AND ENVIRONMENT -----------------------#


#which environmental reading is closest to each krill?
#each row of krill_env corresponds to a krill reading

krill_env <- NULL
for (i in 1:length(krill$datetime)) {
  
  krill_env <- rbind(krill_env, env[which.min(abs(as.numeric(krill$datetime[i] - env$datetime)*24)), ])
  
}

par(mfrow = c(1, 6))
boxplot(krill_env$Sightability ~ whale_present, main = "Sightability")
boxplot(krill_env$SeaState~ whale_present, main = "Sea State")
boxplot(krill_env$Swell~ whale_present, main = "Swell")
boxplot(krill_env$Visibility~ whale_present, main = "Visibility")
boxplot(krill_env$CloudCover~ whale_present, main = "Cloud Cover")
boxplot(krill_env$Intensity~ whale_present, main = "Glare Intensity")

#only significant variables after backwards selection
krill.glm <- glm(whale_number ~ log(arealDen) + krill_env$CloudCover, family = binomial)
summary(krill.glm)

whale_estimate <- as.logical(round(krill.glm$fitted.values))
table(whale_present[!is.na(arealDen)], whale_estimate)

#add environmental variables to hurdle model
krill.hurdle <- hurdle(whale_number ~ log(arealDen) + krill_env$CloudCover, dist = "poisson", zero.dist = "binomial", link = "logit")
summary(krill.hurdle)

whale_estimate <- as.logical(round(krill.hurdle$fitted.values))
table(whale_present[!is.na(arealDen)], whale_estimate)


#how many krill bins within 15km (on same trackline?)

krill_time_difference <- timeDifference(krill, krill)
krill_distance <- apply(krill, 1, FUN = distFromPoint, krill = krill, gps = gps, truePosition=TRUE)

krill_distance[krill_time_difference > 1] <- NA #remove values from another day

#------------------------------- RASTER MAP OF KRILL AND WHALES --------------------------------- #

#for each lat/long on the map, give an effort score
#effort = 0 if no krill cells within horizon distance (13.46km)

#find box size distance (km)
km <- 10
box_x <- round(gcdHF(deg2rad(-67.7), deg2rad(162), deg2rad(-67.7), deg2rad(165.2))/km)
box_y <- round(gcdHF(deg2rad(-67.7), deg2rad(162), deg2rad(-66), deg2rad(162))/km)

location_grid <- raster(ncol = box_x, nrow = box_y, xmn = 162, xmx = 165.2, ymn = -67.7, ymx = -66)

krill$arealDen[krill$arealDen == 0] <- NA
krill_raster <- rasterize(cbind(krill$Longitude, krill$Latitude), location_grid, log(krill$arealDen), fun = mean)
whale_raster <- rasterize(rev(true_lat_long), location_grid, rep(1, nrow(true_lat_long)), fun = sum)
effort       <- rasterize(cbind(gps$Longitude, gps$Latitude), location_grid, gps$bin_time, fun = sum)
sea_state    <- rasterize(cbind(krill$Longitude, krill$Latitude), location_grid, krill_env$SeaState, fun = mean)
sightability <- rasterize(cbind(krill$Longitude, krill$Latitude), location_grid, krill_env$Sightability, fun = mean)
cloud        <- rasterize(cbind(krill$Longitude, krill$Latitude), location_grid, krill_env$CloudCover, fun = mean)

#create raster stack of predictors
predictors <- stack(krill_raster, sea_state, sightability, cloud)
names(predictors) <- c('krill', 'sea_state', 'sightability', 'cloud') 
plot(predictors)

whale_zeros <- getValues(whale_raster)
whale_zeros[!is.na(getValues(krill_raster)) & is.na(getValues(whale_raster))] <- 0
whale_raster <- setValues(whale_raster, whale_zeros)

par(mfrow = c(1, 2))
plot(whale_raster, main = "Total whale sightings")
plot(effort, main = "Effort - time spent in cell (mins)")

#dist from top in km
x <- coordinates(krill_raster)[, 1]
y <- coordinates(krill_raster)[, 2]

cell_x <- colFromX(effort, x)
cell_y <- rowFromY(effort, y)

d <- data.frame(cbind(getValues(whale_raster)/getValues(effort)*60, getValues(predictors), x, y))
d <- data.frame(cbind(d[, 1], apply(d[, c(2:7)], 2, FUN = scale, scale = FALSE)))
names(d) <- c("whales_per_hour", "krill", "sea_state", "sightability", "cloud", "long", "lat")
d$whales_per_hour[!is.na(d$krill) & is.na(d$whales_per_hour)] <- 0
d <- na.omit(d)

whale_pa <- rep(0, nrow(d))
whale_pa[d$whales_per_hour > 0] <- 1


calcPA <- function (model, reference, data) {
  
  #function to print observed vs fitted presence/absence and sensitivity/specificity
  #model: model object
  #reference: reference whale_pa data set
  
  if (class(model)[1] == "hurdle") {
    
    ilogit <- function(x) 1/(1+exp(-x))
    presence_absence <- round(ilogit(data$krill * model$coefficients$zero[1]))
    
    
  } else {
    
    presence_absence <- round(fitted(model))
    presence_absence[presence_absence > 0] <- 1
    
  }
  
  contingency_table <- table(reference, presence_absence, dnn = c("reference", "fitted"))
  
  print(contingency_table)
  
  sens <- round(sensitivity(as.factor(presence_absence), as.factor(reference)), 2)
  spec <- round(specificity(as.factor(presence_absence), as.factor(reference)), 2)
  
  cat("\n", paste("Sensitivity =", sens))
  cat("\n", paste("Specificity =", spec))
  
  return(list("Sensitivity" = sens, "Specificity" = spec))
  
}

#---------------------------------- MODELS ---------------------------------#

#poisson hurdle model
raster.hurdle <- hurdle(round(whales_per_hour) ~ sightability + lat*long + cloud | 
                     krill - 1, dist = "poisson", zero.dist = "binomial", link = "logit", data = d)
summary(raster.hurdle)

hurdle.ss <- calcPA(raster.hurdle, whale_pa, d)

#zero inflated poisson model
raster.zeroinfl <- zeroinfl(round(whales_per_hour)  ~ sightability + sea_state + cloud + krill + lat*long | 
                              krill - 1, dist = "poisson", link = "logit", data = d)
summary(raster.zeroinfl)

zeroinfl.ss <- calcPA(raster.zeroinfl, whale_pa, d)

#poisson glm


raster.glm <- glm(round(whales_per_hour) ~ krill*long + lat*long + sightability + sea_state - 1, family = "poisson", data = d)
summary(raster.glm)

vif(raster.glm)

dispersiontest(raster.glm, alternative ="greater") #test for overdispersion of poisson glm

raster.null <- glm(round(whales_per_hour) ~ 1, family = "poisson", data = d)
anova(raster.glm, raster.null, test = "Chi") #analysis of deviance against the null model

glm.ss <- calcPA(raster.glm, whale_pa, d)

#raster plot of observed vs predicted on the same colour scale

predicted <- effort
predicted_full <- rep(NA, length(getValues(whale_raster)))
predicted_full[as.numeric(names(fitted(raster.glm)))] <- round(fitted(raster.glm))
predicted_full <- predicted_full*getValues(effort)/60
predicted <- setValues(predicted, predicted_full)

par(mfrow = c(1, 2))
plot(whale_raster, col=rev(terrain.colors(ceiling(maxValue(predicted)))), breaks = seq(0, ceiling(maxValue(predicted))))
plot(predicted, col=rev(terrain.colors(ceiling(maxValue(predicted)))), breaks = seq(0, ceiling(maxValue(predicted))))


#negative binomial model

raster.nb <- glm.nb(round(whales_per_hour) ~ krill*lat + lat*long, data = d, maxit = 1000)
summary(raster.nb)

nb.ss <- calcPA(raster.nb, whale_pa, d)

#quasipoisson

raster.qpois <- glm(round(whales_per_hour) ~ krill*long + lat + sightability + sea_state - 1, family = "quasipoisson", data = d)
summary(raster.qpois)

qpois.ss <- calcPA(raster.qpois, whale_pa, d)

#random forest 

raster.rf <- randomForest(round(whales_per_hour) ~ krill + lat + long + sea_state + sightability + cloud, data = d)
raster.rf$importance


#------------------------------ COMPARE MODELS ---------------------------------#


#compare log likelihoods

c("hurdle" = logLik(raster.hurdle), "ZIP" =  logLik(raster.zeroinfl), "Pois" = logLik(raster.glm), "Quasi-Pois" = logLik(raster.qpois), "NB" = logLik(raster.nb))


#compare zero counts

round(c("Obs" = sum(d$whales_per_hour < 1), "hurdle" = sum(predict(raster.hurdle, type = "prob")[, 1]), 
        "ZIP" = sum(predict(raster.zeroinfl, type = "prob")[, 1]),
          "Pois" = sum(dpois(0, fitted(raster.glm))), "Quasi-Pois" = sum(dpois(0, fitted(raster.qpois))), 
        "NB" = sum(dnbinom(0, mu = fitted(raster.nb), size = raster.nb$theta)), 
      "Random Forest" = sum(raster.rf$predicted < 0.5)))


#Observed vs fitted plots

par(mfrow = c(2, 3))

plot(d$whales_per_hour, fitted(raster.hurdle), pch = 19, ylim = c(0, 65), main = "Hurdle")
points(c(0, 100), c(0, 100), col = "red", type = "l")

plot(d$whales_per_hour, fitted(raster.zeroinfl), ylim = c(0, 65), pch = 19, main = "ZIP")
points(c(0, 100), c(0, 100), col = "red", type = "l")

plot(d$whales_per_hour, fitted(raster.glm), pch = 19, main = "Poisson")
points(c(0, 100), c(0, 100), col = "red", type = "l")

plot(d$whales_per_hour, fitted(raster.nb), pch = 19, ylim = c(0, 65), main = "NB")
points(c(0, 100), c(0, 100), col = "red", type = "l")

plot(d$whales_per_hour, fitted(raster.qpois), pch = 19, main = "Quasi-Poisson", ylim = c(0, 65))
points(c(0, 100), c(0, 100), col = "red", type = "l")

plot(d$whales_per_hour, raster.rf$predicted, pch = 19, main = "Random Forest", ylim = c(0, 65))
points(c(0, 100), c(0, 100), col = "red", type = "l")


#plot sensitivity and specificity for each model

ss_type <- data.frame(factor(), numeric(), factor()) 
new_rows <- cbind("model" = rep(c("hurdle", "zeroinfl", "glm", "nb", "qpois"), each = 2), 
                  "value" = rep(c("sensitivity", "specificity"), 5),
                 "sens" = as.numeric(c(unlist(hurdle.ss), unlist(zeroinfl.ss), unlist(glm.ss), unlist(nb.ss), unlist(qpois.ss))))
ss_type <- rbind(ss_type, new_rows)

ggplot(ss_type, aes(value, sens, group = model, size = 2, colour = model)) + 
  theme_bw() + 
  geom_point() + 
  geom_line() + 
  xlab("") + 
  ylab("") + 
  guides(size = FALSE)



# ------------------------- GEOGRAPHICALLY WEIGHTED REGRESSION --------------------------#


d <- data.frame(cbind(getValues(whale_raster), getValues(effort), getValues(predictors), x, y))
d <- data.frame(cbind(d[, c(1, 7, 8)], apply(d[, c(2:6)], 2, FUN = scale, scale = FALSE)))
names(d) <- c("whales", "long", "lat", "effort", "krill", "sea_state", "sightability", "cloud")
d <- na.omit(d)

whale_pa <- rep(0, nrow(d))
whale_pa[d$whales > 0] <- 1

raster.glm <- glm(whales ~ krill + sea_state + cloud, family = "poisson", data = d)
summary(raster.glm)


dispersiontest(raster.glm, alternative ="greater") #test for overdispersion of poisson glm

raster.null <- glm(whales ~ 1, family = "poisson", data = d)
anova(raster.glm, raster.null, test = "Chi") #analysis of deviance against the null model

glm.ss <- calcPA(raster.glm, whale_pa, d)

plot(d$whales, fitted(raster.glm), pch = 19, main = "Poisson")
points(c(0, 100), c(0, 100), col = "red", type = "l")

#test for spatial autocorrelation using Moran's I

dists <- as.matrix(dist(cbind(d$long, d$lat)))
dists.inv <- 1/dists
diag(dists.inv) <- 0

Moran.I(residuals(raster.glm), dists.inv)




dists <- gw.dist(cbind(d$long, d$lat), longlat = TRUE)

sp.data <- SpatialPointsDataFrame(coords <- cbind(d$long, d$lat), data = d)

raster.gwr.bw <- bw.ggwr(whales ~ krill, data = sp.data, adaptive = FALSE, family = "poisson",
                      longlat = TRUE, dMat = dists)

raster.gwr <- gwr.generalised(whales ~ krill, data = sp.data, bw = raster.gwr.bw, adaptive = FALSE, family = "poisson",
          longlat = TRUE)
raster.gwr

gwr.model.fitted <- exp(d$krill*raster.gwr$SDF$krill + raster.gwr$SDF$Intercept)


par(mfrow = c(1, 2))

plot(d$whales, raster.gwr$glm.res$fitted.values, pch = 19, main = "Poisson", ylim = c(0, 8))
points(c(0, 100), c(0, 100), col = "red", type = "l")

plot(d$whales, gwr.model.fitted, pch = 19, main = "Poisson", ylim = c(0, 8))
points(c(0, 100), c(0, 100), col = "red", type = "l")




results <- as.data.frame(raster.gwr$SDF)

ggplot(results, aes(x=coords.x1,y=coords.x2))+geom_point(aes(colour=krill, size = 2)) +
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar", guide_legend(title="Krill")) + 
  guides(size = FALSE) + 
  theme_bw()




