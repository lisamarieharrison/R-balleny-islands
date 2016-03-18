#balleny islands density surface modelling
#author: Lisa-Marie Harrison
#date: 16/03/2016

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/csv")
  source_location <- "C:/Users/43439535/Documents/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Balleny Islands/csv")
  source_location <- "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
}

gps      <- read.csv("GpsData.csv", header = T)
sighting <- read.csv("Sighting.csv", header = T)
env      <- read.csv("Environment.csv", header = T)
effort   <- read.csv("Effort.csv", header = T)
krill    <- read.csv("Krill.csv", header = T)
reticle  <- read.csv("reticle.csv", header = T)
library(chron)
library(ggplot2)
library(geosphere) #destPoint
library(maptools) #gcDestination
library(mapdata)
library(maps)
library(raster) 
library(dsm)
library(mrds)
library(dsm) #density surface model
library(Distance)
library(sp)
library(rgdal)

#source required functions
function_list <- c("gcdHF.R",
                   "deg2rad.R",
                   "rocCurve.R",
                   "draw_map_scale.R",
                   "getFittedGWR.R",
                   "calcPA.R",
                   "distToCell.R",
                   "interpolateWithBarriers.R",
                   "withinTimes.R",
                   "onEffort.R",
                   "sightingDistance.R",
                   "sightingAngle.R",
                   "sightingLatLong.R",
                   "distFromKrill.R",
                   "krillBinTimeDiff.R",
                   "krillWeightedAverage.R",
                   "countIndividalsAroundKrill.R",
                   "removeRasterOverlap.R",
                   "dfFromRaster.R"
)

for (f in function_list) {
  source(paste(source_location, f, sep = ""))
}

#subset sightings to only HB whales (HB = Species 07) & to only MI platform
sighting <- subset(sighting, Species == 07 & Platform == "MI")

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


#------------------------------------ TRUE SIGHTING LOCATION ---------------------------------------#

#find true lat and long of sighting using reticle
#average human eye height on monkey island height is 15.08m

sighting$distance <- unlist(apply(sighting, 1, sightingDistance, reticle = reticle))


#sighting location given specified distance

sighting$angle_true <- apply(sighting, 1, sightingAngle, gps = gps)
true_lat_long <- data.frame(t(apply(sighting, 1, sightingLatLong, gps = gps)))

true_lat_long <- SpatialPoints(na.omit(rev(true_lat_long)), proj4string = CRS("+proj=longlat +datum=WGS84"))
true_lat_long_utm <- spTransform(true_lat_long, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))


#---------------------------- ALLOCATE POINTS TO TRANSECTS ----------------------------------#

direction <- gps$Heading
x <- gps$Longitude
y <- gps$Latitude

plot(krill$Longitude, krill$Latitude, col = "white")
text(krill$Longitude, krill$Latitude, c(1:nrow(krill)), cex = 0.5)

krill$transect <- rep(NA, nrow(krill))
krill$transect[1:55]    <- 1
krill$transect[56:108]  <- 2
krill$transect[109:160] <- 3
krill$transect[161:293] <- 4


# --------------------------- CALCULATE PERPENDICULAR DISTANCES -----------------------------#

#distance to closest point on transect (krill cell)
sighting$angle_true <- apply(sighting, 1, sightingAngle, gps = gps)
true_lat_long <- data.frame(t(apply(sighting, 1, sightingLatLong, gps = gps)))
distance   <- apply(true_lat_long, 1, FUN = distFromKrill, krill = krill, gps = gps, truePosition = TRUE)

closest_bin <- apply(distance, 2, which.min)
closest_bin[is.na((closest_bin == Inf))] <- NA
closest_bin <- unlist(closest_bin)
closest_distance <- apply(distance, 2, min)*1000

#remove error values
krill$arealDen[krill$arealDen == 0] <- 0.0001

obs_count <- rep(0, nrow(krill))
obs_count[as.numeric(names(table(closest_bin)))] <- table(closest_bin)


#------------------------------ DENSITY SURFACE MODEL ----------------------------------------#


balleny_map <- map("world2Hires", regions=c("Antarctica:Young Island", "Antarctica:Buckle Island", "Antarctica:Sturge Island"), plot = FALSE)
balleny_poly <- map2SpatialPolygons(balleny_map, IDs = balleny_map$names, proj4string=CRS("+proj=longlat +datum=WGS84"))
balleny_poly_utm <- spTransform(balleny_poly, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))


#fit detection function

distdata <- data.frame(cbind(c(1:nrow(sighting)), sighting$BestNumber, sighting$distance*1000, rep(1, nrow(sighting)), gps$Latitude[match(sighting$GpsIndex, gps$Index)], gps$Longitude[match(sighting$GpsIndex, gps$Index)]))
colnames(distdata) <- c("object", "size", "distance", "detected", "latitude", "longitude")

det_function <- ds(distdata, max(distdata$distance), key="hr", adjustment=NULL)
#det_function_size <-ds(distdata, max(distdata$distance), formula=~as.factor(size), key="hr", adjustment=NULL)
summary(det_function)
plot(det_function)


xy <- SpatialPoints(cbind(krill$Longitude, krill$Latitude))
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
res <- spTransform(xy, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))



segdata <- data.frame(cbind(krill$Longitude, krill$Latitude, coordinates(res), krill$Distance_vl, krill$transect, c(1:nrow(krill)), log(krill$arealDen), obs_count))
colnames(segdata) <- c("longitude", "latitude", "x", "y", "Effort", "Transect.Label", "Sample.Label", "krill", "number")


obsdata <- data.frame(cbind(c(1:nrow(sighting)), closest_bin, sighting$BestNumber, sighting$distance*1000))
names(obsdata) <- c("object", "Sample.Label", "size", "distance")
obsdata <- na.omit(obsdata)


whale.dsm <- dsm(formula = count ~ s(x, y) + krill, ddf.obj = det_function, family = "poisson", segment.data = segdata, observation.data = obsdata, method="REML")
summary(whale.dsm)

#plot relative counts over the smooth space
vis.gam(whale.dsm, plot.type="contour", view = c("x","y"), too.far = 0.06, asp = 1, type = "response", contour.col = "black", n.grid = 100)
plot(balleny_poly_utm, add = TRUE, col = "grey")
points(true_lat_long_utm, col = "blue", pch = 19)

#plot count relationship to krill

#plot observed vs fitted
plot(na.omit(segdata)$number, whale.dsm$fitted.values)







