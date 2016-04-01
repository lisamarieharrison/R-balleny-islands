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
library(plyr) #join
library(AER) #dispersiontest
library(rgeos) #gArea

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
                   "dfFromRaster.R",
                   "krillToGrid.R",
                   "grid_plot_obj.R",
                   "check_cols.R"
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

#effort where MI observers pdat_loc_utment
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

#which environmental reading is closest to each krill?
#each row of krill_env cordat_loc_utmponds to a krill reading
krill_env <- NULL
for (i in 1:length(krill$datetime)) {
  
  krill_env <- rbind(krill_env, env[which.min(abs(as.numeric(krill$datetime[i] - env$datetime)*24)), ])
  
}

#------------------------------------ TRUE SIGHTING LOCATION ---------------------------------------#

#find true lat and long of sighting using reticle
#average human eye height on monkey island height is 15.08m

sighting$distance <- unlist(apply(sighting, 1, sightingDistance, reticle = reticle))

#remove sightings > 10km away because no reticle between 6.5km - 13.8km so distance inacurate
sighting <- sighting[sighting$distance < 10, ]

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
time_difference <- krillBinTimeDiff(sighting, krill)
distance   <- apply(true_lat_long, 1, FUN = distFromKrill, krill = krill, gps = gps, truePosition = TRUE)
distance[time_difference > 1] <- NA # remove distances over 1 hour away from sighting

closest_bin <- apply(distance, 2, which.min)
closest_bin[is.na((closest_bin == Inf))] <- NA
closest_bin <- unlist(closest_bin)
closest_distance <- apply(distance, 2, min)*1000

#remove error values
krill$arealDen[krill$arealDen == 0] <- 0.0001

obs_count <- rep(0, nrow(krill))
obs_count[as.numeric(names(table(closest_bin)))] <- table(closest_bin)


# ----------------------------- DENSITY SURFACE MODEL -------------------------------- #

#read csv file of coordinates for each island
for (i in 1:3) {
  
  b <- read.csv(paste("C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/polygons/Balleny_", i, ".csv", sep = ""), header = F)
  x <- unlist(b[seq(1, length(b), by = 3)])
  y <- unlist(b[seq(2, length(b), by = 3)])
  assign(paste("b", i, "_poly", sep = ""), Polygon(cbind(x, y), hole = FALSE))
  
}

balleny_poly     <- SpatialPolygons(list(Polygons(list(b1_poly), ID = 1), Polygons(list(b2_poly), ID = 2), Polygons(list(b3_poly), ID = 3)), proj4string = CRS("+proj=longlat +datum=WGS84"))
balleny_poly_utm <- spTransform(balleny_poly, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))
balleny_ggplot   <- fortify(balleny_poly_utm, region="id") #df for ggplot

#convert data locations to spatial points
dat_loc     <- SpatialPoints(cbind(krill$Longitude, krill$Latitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
dat_loc_utm <- spTransform(dat_loc, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))


#fit detection function

segdata <- data.frame(cbind(krill$Longitude, krill$Latitude, coordinates(dat_loc_utm), krill$Distance_vl, krill$transect, c(1:nrow(krill)), log(krill$arealDen), obs_count, krill_env$CloudCover, krill_env$SeaState, krill_env$Sightability))
colnames(segdata) <- c("longitude", "latitude", "x", "y", "Effort", "Transect.Label", "Sample.Label", "krill", "number", "cloud", "sea_state", "sightability")

obsdata <- data.frame(cbind(c(1:nrow(sighting)), closest_bin, segdata$Transect.Label[closest_bin], sighting$BestNumber, sighting$distance*1000))
names(obsdata) <- c("object", "Sample.Label", "Transect.Label", "size", "distance")

distdata <- data.frame(cbind(c(1:nrow(sighting)), sighting$BestNumber, sighting$distance*1000, rep(1, nrow(sighting)), 
            gps$Latitude[match(sighting$GpsIndex, gps$Index)], gps$Longitude[match(sighting$GpsIndex, gps$Index)], 
            obsdata$Sample.Label, obsdata$Transect.Label, krill_env$SeaState[closest_bin], krill_env$SeaState[closest_bin]))
colnames(distdata) <- c("object", "size", "distance", "detected", "latitude", "longitude", "Sample.Label", "Transect.Label", "sea_state", "cloud")

obsdata  <- na.omit(obsdata)
distdata <- na.omit(distdata)
segdata  <- na.omit(segdata)

#for ds function
region.table <- segdata[6]
region.table$Area <- segdata$Effort*(7000 - 200)*2
region.table <- aggregate(region.table$Area, by = list(region.table$Transect.Label), FUN = "sum")
names(region.table) <- c("Region.Label", "Area")

sample.table <- segdata[5:7]
names(sample.table) <- c("Effort", "Region.Label", "Sample.Label")

obs.table <- obsdata[1:3]
names(obs.table) <- c("object", "Sample.Label", "Region.Label")

#using ds
det_function <- ds(data = distdata, truncation = list(left = 200, right = 7000) , key="hn", adjustment=NULL, sample.table = sample.table, region.table = region.table, obs.table = obs.table)
det_function_size <- ds(distdata, truncation = list(left = 200, right = 7000) , formula=~size, key="hn", adjustment=NULL, sample.table = sample.table, region.table = region.table, obs.table = obs.table)
summary(det_function_size)
plot(det_function)

#using mrds
det_function <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
               data = distdata, meta.data = list(left = 200, width = 7000))
#summary(det_function)
#plot(det_function)

#calculate area of each segment using length of segment
segment.area <- segdata$Effort*(7000 - 200)*2

whale.dsm <- dsm(formula = D ~ s(x, y, k = 10) + krill, family = tw(), ddf.obj = det_function, segment.data = segdata, observation.data = obsdata, method="REML", segment.area = segment.area)
summary(whale.dsm)

#plot relative counts over the smooth space
vis.gam(whale.dsm, plot.type="contour", view = c("x","y"), too.far = 0.06, asp = 1, type = "dat_loc_utmponse", contour.col = "black", n.grid = 100)
plot(balleny_poly_utm, add = TRUE, col = "grey")
points(true_lat_long_utm, col = "blue", pch = 19)

#goodness of fit
gam.check(whale.dsm)

#check overdispersion if poisson glm
dispersiontest(whale.dsm, alternative ="greater")

#check spatial autocorrelation
dsm.cor(whale.dsm, max.lag = 10, Segment.Label="Sample.Label")

#AIC check
pchisq(493.1144 -  491.9856, 1, lower.tail=FALSE)


# ------------------------------ SURVEY AREA POLYGON -------------------------------- #

#calculate convex hull around points
ch <- chull(cbind(segdata$x, segdata$y))
coords <- cbind(segdata$x, segdata$y)[c(ch, ch[1]), ]  # closed polygon

pred.polys <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)), proj4string = CRS("+proj=utm +zone=58 +south +ellps=WGS84"))

grid <- raster(extent(pred.polys))

# Choose its res (m)
res(grid) <- 10000

# Make the grid have the same coordinate reference system (CRS) as the shapefile.
proj4string(grid) <- proj4string(pred.polys)

#get percentage of cells overlapped by islands
overlap_poly <- getValues(rasterize(balleny_poly_utm, grid, getCover = TRUE))
grid <- setValues(grid, overlap_poly)

# Transform this raster into a polygon to create grid
gridpolygon <- rasterToPolygons(grid)

#Intersect with survey area
survey.grid <- intersect(pred.polys, gridpolygon)

survey_area <- gArea(pred.polys) - gArea(balleny_poly_utm) #area in m2 minus island area

# ------------------------------ SOAP FILM SMOOTHER ---------------------------- #


#soap smoother to remove island
island.hole <- gDifference(survey.grid, balleny_poly_utm)
island.grid <- intersect(island.hole, gridpolygon)
knot_points <- list(x = coordinates(island.grid)[, 1], y= coordinates(island.grid)[, 2])
soap.knots  <- make.soapgrid(knot_points, c(10, 10))

#increase survey area by 10km

grid <- raster(extent(gBuffer(survey.grid, width = 10000)))
# Choose its dat_loc_utmolution (m)
res(grid) <- 10000

# Make the grid have the same coordinate reference system (CRS) as the shapefile.
proj4string(grid)<-proj4string(gBuffer(survey.grid, width = 10000))

#get percentage of cells overlapped by islands
overlap_poly <- getValues(rasterize(balleny_poly_utm, grid, getCover = TRUE))
grid <- setValues(grid, overlap_poly)

# Transform this raster into a polygon to create grid
gridpolygon <- rasterToPolygons(grid)

survey.grid.large <- intersect(gBuffer(survey.grid, width = 10000), gridpolygon)

#remove boundary points
ch     <- chull(coordinates(survey.grid.large))
coords <- coordinates(survey.grid.large)[c(ch, ch[1]), ] 

#calculate area of each cell (m)
grid_cell_area <- rep((res(grid)[1])^2, nrow(coordinates(survey.grid)))*(1-survey.grid$layer/100)

#calculate weighted krill around each point
krill_mean <- apply(coordinates(survey.grid), 1, krillToGrid, threshold = res(grid)[1]/1000, krill_mat = segdata)

#bnd is list of islands boundaries (survey area and 3 islands) which can't overlap
bnd <- list(xy.coords(coords), xy.coords(fortify(balleny_poly_utm[1])[, 1:2]), xy.coords(fortify(balleny_poly_utm[2])[, 1:2]), xy.coords(fortify(balleny_poly_utm[3])[, 1:2]))

#remove knots inside islands
x <- soap.knots[, 1]
y <- soap.knots[, 2]
soap.knots <- soap.knots[inSide(bnd, x, y), ]

#check data format is correct
check.cols(ddf.obj = det_function, segment.data = segdata, observation.data = obsdata, segment.area = segment.area)

whale.dsm <- dsm(D ~ s(x, y, bs="so", k = 5, xt=list(bnd=bnd)) + krill, family = tw(), ddf.obj = det_function, 
                 segment.data = segdata, observation.data = obsdata, method="REML", segment.area = segment.area, 
                 knots = soap.knots)
summary(whale.dsm)


#plot relative counts over the smooth space
vis.gam(whale.dsm, plot.type="contour", view = c("x","y"), too.far = 0.1, asp = 1, type = "response", contour.col = "black", n.grid = 100)
plot(balleny_poly_utm, add = TRUE, col = "grey")
points(true_lat_long_utm, pch = 19)

#goodness of fit
gam.check(whale.dsm)
rqgam.check(whale.dsm) #randomised quantile residuals

#check spatial autocorrelation
dsm.cor(whale.dsm, max.lag = 10, Segment.Label="Sample.Label")


# ---------------------------- ABUNDANCE ESTIMATION ---------------------------#

#data frame of prediction locations
preddata <- data.frame(cbind(coordinates(survey.grid), grid_cell_area, krill_mean, res(grid)[1]), rep(c(1:12), each = 8))
colnames(preddata) <- c("x", "y", "area", "krill", "Effort")

#calculate predicted values
if (all.vars(whale.dsm$formula)[1] == "D") {
  
  #if density model, predict densities in whales/km^2
  whale_pred <- c(predict(whale.dsm, preddata, off.set = 0))/1e-6
  total_individuals <- mean(na.omit(whale_pred))*survey_area*1e-6
  cv <- c(predict(whale.dsm, preddata, off.set = 0, se.fit = TRUE)$se.fit)/(whale_pred*1e-6)
  ddf.cv <- summary(whale.dsm$ddf)$average.p.se/summary(whale.dsm$ddf)$average.p

} else {
  
  #if abundance model, predict abundance in each cell
  whale_pred <- c(predict(whale.dsm, preddata, off.set = preddata$area))
  total_individuals <- sum(na.omit(whale_pred))
  
}

p <- ggplot() + grid_plot_obj(fill = whale_pred, name = all.vars(whale.dsm$formula)[1], sp = survey.grid) +
  geom_polygon(data=balleny_ggplot, aes(x=long, y=lat, group=id), color="black", fill = "grey")
p

p <- ggplot() + grid_plot_obj(fill = cv + ddf.cv, name = "CV", sp = survey.grid) +
  geom_polygon(data=balleny_ggplot, aes(x=long, y=lat, group=id), color="black", fill = "grey") + 
  geom_point(aes(x = Longitude, y = Latitude), data = data.frame(coordinates(true_lat_long_utm)))
p


# -------------------------- VARIANCE ESTIMATION ---------------------------#

#dsm.var.prop can't handle NA values in preddata so need to na.omit and keep track of data location with prediction.points
preddata_na <- na.omit(preddata)
prediction_points <- which(rowSums(is.na(preddata)) == 0)
preddata.varprop <- split(preddata_na, 1:nrow(preddata_na))

if (all.vars(whale.dsm$formula)[1] == "D") {
  
  #if density model, predict densities in whales/km^2
  dsm.xy.varprop <- dsm.var.prop(whale.dsm, pred.data = preddata.varprop, off.set = 0)

  dsm.xy.varprop <- dsm.var.gam(whale.dsm, pred.data = preddata_na, off.set = 0)
  

  } else {
  
  #if abundance model, predict abundance in each cell
  dsm.xy.varprop <- dsm.var.prop(whale.dsm, pred.data = preddata.varprop, off.set = preddata_na$area)
  
  dsm.xy.varprop <- dsm.var.gam(whale.dsm, pred.data = preddata_na, off.set = preddata_na$area)
}


pred_var <- rep(NA, nrow(preddata))
pred_var[prediction_points] <- dsm.xy.varprop$pred.var

pred <- rep(NA, nrow(preddata))
pred[prediction_points] <- dsm.xy.varprop$pred

p <- ggplot() + grid_plot_obj(sqrt(pred_var)/unlist(pred), "CV", sp = survey.grid) +
  geom_polygon(data=balleny_ggplot, aes(x=long, y=lat, group=id), color="black", fill = "grey")
p






