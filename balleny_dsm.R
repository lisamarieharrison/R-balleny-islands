#balleny islands density surface modelling
#author: Lisa-Marie Harrison
#date: 16/03/2016

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "~/Lisa/phd/Balleny Islands/csv")
  source_location <- "~/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
} else {
  setwd(dir = "~/phd/southern ocean/Balleny Islands/csv")
  source_location <- "~/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
}

gps      <- read.csv("GpsData.csv", header = T)
sighting <- read.csv("Sighting.csv", header = T)
env      <- read.csv("Environment.csv", header = T)
effort   <- read.csv("Effort.csv", header = T)
krill    <- read.csv("CombinedKrillDen.csv", header = T)
reticle  <- read.csv("reticle.csv", header = T)
under    <- read.csv("das-data_2015-01-24_2015-03-12.csv", header = T)

library(chron)
library(ggplot2)
library(geosphere) #destPoint
library(maptools) #gcDestination
library(mapdata)
library(maps)
library(raster) 
library(mrds)
library(dsm) #density surface model
library(Distance)
library(sp)
library(plyr) #join
library(AER) #dispersiontest
library(rgeos) #gArea
library(raadsync)
library(raadtools)
library(rgdal) #process hdf files
library(rhdf5) #read hdf files

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
                   "envToGrid.R",
                   "grid_plot_obj.R",
                   "check_cols.R",
                   "reticleDistances.R",
                   "grid_plot_obj_NA.R",
                   "multiplot.R",
                   "draw_map_scale"
)

for (f in function_list) {
  source(paste(source_location, f, sep = ""))
}

#subset sightings to only HB whales (HB = Species 07) & to only MI platform
sighting <- subset(sighting, Species == 07 & Platform == "MI")

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


#which environmental reading is closest to each krill?
#each row of krill_env coresponds to a krill reading
krill_env <- NULL
for (i in 1:length(krill$datetime)) {
  
  krill_env <- rbind(krill_env, env[which.min(abs(as.numeric(krill$datetime[i] - env$datetime)*24)), ])
  
}


# ------------------------------ TOTAL ON EFFORT TIME ----------------------------------- #

#on-effort time (hr)
sum(diff(effort$datetime)[seq(1, nrow(effort), by = 2)])*24


# --------------------------------- SURVEY TRACKLINE ------------------------------------ #

gps_lat_long <- SpatialPoints(na.omit(gps[, c("Longitude", "Latitude")]), proj4string = CRS("+proj=longlat +datum=WGS84"))
gps_lat_long_utm <- spTransform(gps_lat_long, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))


plot(balleny_poly, col = "grey", xlim = c(162, 165.5), ylim = c(-67.7, -66), xlab = "Longitude", ylab = "Latitude")
points(gps_lat_long, pch = 19)

axis(1, seq(162, 165.5, length.out = 5), round(seq(162, 165.5, length.out = 5), 2))

axis(2, seq(-67.7, -66, length.out = 5), round(seq(-67.7, -66, length.out = 5), 2))


#------------------------------------ TRUE SIGHTING LOCATION ---------------------------------------#

#find true lat and long of sighting using reticle
#average human eye height on monkey island height is 15.08m

sighting$distance <- unlist(apply(sighting, 1, sightingDistance, reticle = reticle))

#remove sightings > 10km away because no reticle between 6.5km - 13.8km so distance inacurate
#sighting <- sighting[sighting$distance < 10, ]

#sighting location given specified distance

sighting$angle_true <- apply(sighting, 1, sightingAngle, gps = gps)
true_lat_long <- data.frame(t(apply(sighting, 1, sightingLatLong, gps = gps)))

true_lat_long <- SpatialPoints(na.omit(rev(true_lat_long)), proj4string = CRS("+proj=longlat +datum=WGS84"))
true_lat_long_utm <- spTransform(true_lat_long, CRS("+proj=utm +zone=58 +south +ellps=WGS84"))


#-------------------------------- CALCULATE DISTANCE BINS -----------------------------------#


reticle_dists <- rev(unlist(apply(reticle, 1, reticleDistances)))

left_bin <- 200 #left truncation distance
for (i in 2:length(reticle_dists)) {
  
  left_bin[i] <- (reticle_dists[i] + (reticle_dists[i - 1] - reticle_dists[i])/2)*1000
  
}
left_bin <- c(left_bin, 13800) #right truncation distance

#---------------------------- ALLOCATE POINTS TO TRANSECTS ----------------------------------#

direction <- gps$Heading
x <- gps$Longitude
y <- gps$Latitude

plot(krill$Longitude, krill$Latitude, col = "white")
text(krill$Longitude, krill$Latitude, c(1:nrow(krill)), cex = 0.5)

krill$transect <- rep(NA, nrow(krill))
krill$transect[1:60]    <- 1
krill$transect[61:122]  <- 2
krill$transect[123:194] <- 3
krill$transect[195:nrow(krill)] <- 4


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

obs_count <- rep(0, nrow(krill))
obs_count[as.numeric(names(table(closest_bin)))] <- table(closest_bin)

#-------------------------------- UNDERWAY DATA -------------------------------------#

under_sp <- SpatialPoints(na.omit(cbind(under$GP500_GPLong, under$GP500_GPLat)), proj4string = CRS("+proj=longlat +datum=WGS84"))
under_sp <- spTransform(under_sp, CRSobj = CRS("+proj=utm +zone=58 +south +ellps=WGS84"))
under <- cbind(under[!is.na(under$GP500_GPLong) & !is.na(under$GP500_GPLat), ], coordinates(under_sp))
colnames(under)[62:63] <- c("x", "y")

under$datetime <- chron(dates. = substr(under$utc, 1, 10), times. = substr(under$utc, 12, 19), format = c(dates. = "y-m-d", times. = "h:m:s"), out.format = c(dates = "d/m/y", times = "h:m:s"))
under      <- subset(under, under$datetime >= min(krill$datetime) & under$datetime <= max(krill$datetime))
under$SB21_SB21sal[under$SB21_SB21sal < 26] <- NA #salinity error values
under$EK60_EK60dbt_38[under$EK60_EK60dbt_38 == 0] <- NA #depth error values
under$SB21_SB21dens[under$SB21_SB21dens < 26] <- NA #density error values

#average underway data within krill bins

under_env <- under[, 5:60]
krill_underway <- data.frame()
for (i in 1:nrow(krill)) {
  
  bin_start <- krill$datetime[i] - 5/3600
  bin_end   <- krill$datetime[i] + 5/3600
  
  krill_underway <- rbind(krill_underway, colMeans(under_env[under$datetime >= bin_start & under$datetime <= bin_end, ]))
  
}
colnames(krill_underway) <- colnames(under[, 5:60])


# ----------------------------- DENSITY SURFACE MODEL -------------------------------- #

#read csv file of coordinates for each island
for (i in 1:3) {
  
  b <- read.csv(paste("~/Lisa/phd/Balleny Islands/polygons/Balleny_", i, ".csv", sep = ""), header = F)
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



#-------------- ATTENUATION FUNCTION FOR EFFORT -----------------------#

front_coords <- matrix(NA, ncol = 2, nrow = nrow(krill))
back_coords <- matrix(NA, ncol = 2, nrow = nrow(krill))

for (i in 1:nrow(krill)) {
  
  angle <- gps$Heading[which.min(abs(gps$datetime -  krill$datetime[i]))]
  
  if (angle < 180) {
    
    angle_back <- angle + 180
    
  } else {
    
    angle_back <- angle - 180
    
  }
  
  distance <- krill$integrationInterval.m[1]/2
  
  coords_fwd <- destPoint(p = c(krill$Longitude[i], krill$Latitude[i]),
                          b = angle, d = distance)
  
  coords_back <- destPoint(p = c(krill$Longitude[i], krill$Latitude[i]),
                           b = angle_back, d = distance)
  
  
  coords_fwd <- SpatialPoints(coords_fwd, CRS(proj4string(balleny_poly)))
  coords_fwd <- spTransform(coords_fwd, CRS(proj4string(balleny_poly_utm)))
  
  coords_back <- SpatialPoints(coords_back, CRS(proj4string(balleny_poly)))
  coords_back <- spTransform(coords_back, CRS(proj4string(balleny_poly_utm)))
  
  front_coords[i, ] <- coordinates(coords_fwd)
  back_coords[i, ]  <- coordinates(coords_back)
  
}

percent <- NULL
for (i in 1:nrow(krill)) {
  
  
  line <-  Line(matrix(c(front_coords[i, ], back_coords[i, ]), nrow = 2, byrow = T))
  
  lines <- Lines(list(line), ID = "l")
  
  spatial_line <- SpatialLines(list(lines), proj4string = CRS(proj4string(balleny_poly_utm)))
  
  large_line <- gBuffer(spatial_line, width = 12000, capStyle = "FLAT")
  
  suppressWarnings(overlap <- intersect(balleny_poly_utm, large_line))
  
  if (is.null(overlap)) {
    overlap_area <- 0
  } else {
    overlap_area <- area(overlap)
  }
  
  percent[i] <- overlap_area/area(large_line)
  
}


#fit detection function

segdata <- data.frame("longitude" = krill$Longitude, "latitude" = krill$Latitude, "x" = coordinates(dat_loc_utm)[, 1], "y" = coordinates(dat_loc_utm)[, 2], 
                      "Effort" = krill$integrationInterval.m, "Transect.Label" = krill$transect, "Sample.Label" = c(1:nrow(krill)), 
                      "krill" = krill$krillArealDen.gm2, "number" = obs_count, "cloud" = krill_env$CloudCover, "sea_state" = krill_env$SeaState, 
                      "sightability" = krill_env$Sightability, "SST" = as.numeric(as.character(krill_env$SST)), "datetime" = krill$datetime,
                      "salinity" = krill_underway$SB21_SB21sal, "bottom_depth" = krill_underway$ES60_ES60dbt, "density" = krill_underway$SB21_SB21dens,
                      "chl" = krill_underway$TRIPLET_TripletChl, "overlap" = percent)

obsdata <- data.frame(cbind(c(1:nrow(sighting)), closest_bin, segdata$Transect.Label[closest_bin], sighting$BestNumber, sighting$distance*1000))
names(obsdata) <- c("object", "Sample.Label", "Transect.Label", "size", "distance")

distdata <- data.frame(cbind(c(1:nrow(sighting)), sighting$BestNumber, sighting$distance*1000, rep(1, nrow(sighting)), 
                             gps$Latitude[match(sighting$GpsIndex, gps$Index)], gps$Longitude[match(sighting$GpsIndex, gps$Index)], 
                             obsdata$Sample.Label, obsdata$Transect.Label, krill_env$SeaState[closest_bin], krill_env$CloudCover[closest_bin], krill_env$Sightability[closest_bin]))
colnames(distdata) <- c("object", "size", "distance", "detected", "latitude", "longitude", "Sample.Label", "Transect.Label", "sea_state", "cloud", "sightability")


obsdata  <- na.omit(obsdata)
distdata <- na.omit(distdata)

#for ds function
region.table <- segdata[6]
region.table$Area <- segdata$Effort*(13800 - 200)*2
region.table <- aggregate(region.table$Area, by = list(region.table$Transect.Label), FUN = "sum")
names(region.table) <- c("Region.Label", "Area")

sample.table <- segdata[5:7]
names(sample.table) <- c("Effort", "Region.Label", "Sample.Label")

obs.table <- obsdata[1:3]
names(obs.table) <- c("object", "Sample.Label", "Region.Label")

#using ds
det_function <- ds(data = distdata, truncation = list(left = 200, right = 13800), cutpoints = left_bin, key="hn", adjustment=NULL)
det_function_size <- ds(distdata, truncation = list(left = 200, right = 13800), cutpoints = left_bin, formula=~size, key="hn", adjustment=NULL, sample.table = sample.table, region.table = region.table, obs.table = obs.table)

# ------------------------------ SURVEY AREA POLYGON -------------------------------- #

#calculate area of each segment using length of segment
segment.area <- segdata$Effort*(13800 - 200)*2*(1-segdata$overlap)

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

#increase survey area by 10km

grid <- raster(extent(gBuffer(survey.grid, width = 10000)))
# Choose its dat_loc_utm (m)
res(grid) <- 10000

# Make the grid have the same coordinate reference system (CRS) as the shapefile.
proj4string(grid)<-proj4string(gBuffer(survey.grid, width = res(grid)[1]))

#get percentage of cells overlapped by islands
overlap_poly <- getValues(rasterize(balleny_poly_utm, grid, getCover = TRUE))
grid <- setValues(grid, overlap_poly)

# Transform this raster into a polygon to create grid
gridpolygon <- rasterToPolygons(grid)

survey.grid.large <- intersect(gBuffer(survey.grid, width = res(grid)[1]), gridpolygon)

# ------------------------------ SOAP FILM SMOOTHER ---------------------------- #


#soap smoother to remove island
island.hole <- gDifference(survey.grid, balleny_poly_utm)
island.grid <- intersect(island.hole, gridpolygon)
knot_points <- list(x = coordinates(island.grid)[, 1], y= coordinates(island.grid)[, 2])
soap.knots  <- make.soapgrid(knot_points, c(11, 11))

#remove boundary points
ch     <- chull(coordinates(survey.grid.large))
coords <- coordinates(survey.grid.large)[c(ch, ch[1]), ] 

#calculate area of each cell (m)
grid_cell_area <- rep((res(grid)[1])^2, nrow(coordinates(survey.grid)))*(1-survey.grid$layer/100)

#calculate weighted krill around each point
krill_mean <- apply(coordinates(survey.grid), 1, envToGrid, threshold = res(grid)[1]/1000, data_frame = segdata, variable = "krill")
cloud_mean <- apply(coordinates(survey.grid), 1, envToGrid, threshold = res(grid)[1]/1000, data_frame = segdata, variable = "cloud")
SST_mean   <- apply(coordinates(survey.grid), 1, envToGrid, threshold = res(grid)[1]/1000, data_frame = segdata, variable = "SST")
salinity_mean   <- apply(coordinates(survey.grid), 1, envToGrid, threshold = res(grid)[1]/1000, data_frame = under, variable = "SB21_SB21sal")
depth_mean      <- apply(coordinates(survey.grid), 1, envToGrid, threshold = res(grid)[1]/1000, data_frame = under, variable = "EK60_EK60dbt_38")
chl_mean      <- apply(coordinates(survey.grid), 1, envToGrid, threshold = res(grid)[1]/1000, data_frame = under, variable = "TRIPLET_TripletChl")


#bnd is list of islands boundaries (survey area and 3 islands) which can't overlap
bnd <- list(xy.coords(coords), xy.coords(fortify(balleny_poly_utm[1])[, 1:2]), xy.coords(fortify(balleny_poly_utm[2])[, 1:2]), xy.coords(fortify(balleny_poly_utm[3])[, 1:2]))

#remove knots inside islands
x <- soap.knots[, 1]
y <- soap.knots[, 2]
soap.knots <- soap.knots[inSide(bnd, x, y), ]

#check data format is correct
check.cols(ddf.obj = det_function, segment.data = segdata, observation.data = obsdata, segment.area = segment.area)

segdata$bottom_depth[segdata$bottom_depth == 0] <- NA

whale.dsm <- dsm(Nhat ~ s(x, y, bs="sw", xt=list(bnd=bnd)) + s(krill) + chl + s(salinity) + s(bottom_depth), ddf.obj = det_function_size, 
                 segment.data = segdata, observation.data = obsdata, method = "REML", segment.area = segment.area, knots = soap.knots)
summary(whale.dsm)


segdata$count <- 0

for (i in 1:nrow(obsdata)) {
  
  seg <- (segdata$Transect.Label == obsdata$Transect.Label[i] & segdata$Sample.Label == obsdata$Sample.Label[i])
  
  segdata$count[seg] <- segdata$count[seg] + 1
  
}


#plot relative counts over the smooth space
vis.gam(whale.dsm, plot.type="contour", view = c("x","y"), too.far = 0.1, asp = 1, type = "response", contour.col = "black", n.grid = 100)
plot(balleny_poly_utm, add = TRUE, col = "grey")
points(true_lat_long_utm, pch = 19)

#check spatial autocorrelation
dsm.cor(whale.dsm, max.lag = 10, Segment.Label="Sample.Label")

#contour plot of x, y with colour scale bar
plot_dat <- plot.gam(whale.dsm, select = 1)
fill_col <- plot_dat[[1]]$fit
x <- plot_dat[[1]]$xscale
y <- plot_dat[[1]]$yscale

image.plot(x, y, fill_col, col = heat.colors(100), xlab = "Easting", ylab = "Northing")
contour(x, y, fill_col, add = TRUE)

# ---------------------------- ABUNDANCE ESTIMATION ---------------------------#

#data frame of prediction locations
preddata <- data.frame(cbind(coordinates(survey.grid), grid_cell_area, krill_mean, res(grid)[1], cloud_mean, SST_mean, salinity_mean, depth_mean, chl_mean, 0, 1))
colnames(preddata) <- c("x", "y", "area", "krill", "Effort", "cloud", "SST", "salinity", "bottom_depth", "chl", "mx0_bottom_depth", "idx0_bottom_depth")

#calculate predicted values

#if density model, predict densities in whales/km^2

if (class(whale.dsm)[2] == "gamm") {
  ddf.obj   <- whale.dsm$gam$ddf
  model.obj <- whale.dsm$gam
} else {
  ddf.obj   <- whale.dsm$ddf
  model.obj <- whale.dsm
}

if (all.vars(model.obj$formula)[1] == "D") {
  whale_pred <- c(predict(whale.dsm, preddata, off.set = 0, se.fit = TRUE)$fit/1e-6)
  print(total_individuals <- mean(na.omit(whale_pred))*survey_area*1e-6)
  cv <- c(predict(whale.dsm, preddata, off.set = 0, se.fit = TRUE)$se.fit/1e-6/whale_pred)
} else {
  whale_pred <- c(predict(whale.dsm, preddata, off.set = preddata$area))
  print(total_individuals <- sum(na.omit(whale_pred)))
  cv <- c(predict(whale.dsm, preddata, off.set = preddata$area, se.fit = TRUE)$se.fit)/whale_pred
}
ddf.cv <- summary(model.obj$ddf)$average.p.se/summary(model.obj$ddf)$average.p

p1 <- ggplot() + grid_plot_obj_NA(fill = whale_pred, name = "Estimate", sp = survey.grid) +
  geom_polygon(data=balleny_ggplot, aes(x=long, y=lat, group=id), color="black", fill = "grey") +
  theme_bw() +
  theme(text =  element_text(size = 25)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = 350000) + 
  geom_hline(yintercept = 2490000) +
  xlab("Easting") +
  ylab("Northing") +
  annotate("text", x = 495000, y = 2655000, label = "(a)", size = 10)

p2 <- ggplot() + grid_plot_obj_NA(fill = cv + ddf.cv, name = "CV", sp = survey.grid) +
  geom_polygon(data=balleny_ggplot, aes(x=long, y=lat, group=id), color="black", fill = "grey") + 
  theme_bw() +
  theme(text =  element_text(size = 25)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  geom_vline(xintercept = 350000) + 
  geom_hline(yintercept = 2490000) +
  xlab("Easting") +
  ylab("") +
  annotate("text", x = 495000, y = 2655000, label = "(b)", size = 10)

multiplot(p1, p2, cols = 2)

#SE of prediction
sqrt(sum(na.omit(cv*whale_pred)^2))

#------------------------------- SEA ICE DATA AAD --------------------------------#

#get sea ice percentage coverage from AAD data centre using raadtools

cfg <- read_repo_config(system.file("extdata", "raad_repo_config.json", package= "raadsync"))
ice_index <- 1
cfg$do_sync <- seq(nrow(cfg)) == ice_index

## limit the data to only Feb 2015
cfg$method_flags[1] <- paste0(cfg$method_flags[1], " --accept=\"*nt_201502*\"")

## specify local repository location
my_datadir <- normalizePath("~/Lisa/phd/Balleny Islands/remote data/sea ice", "/")
options(default.datadir = my_datadir)
cfg$local_file_root <- file.path(my_datadir, "data")

for (i in 3:6) {
  
  ice <- readice(date=chron(dates. = paste0("2015/02/0", i), format = "Y/m/d"))
  
  balleny_ice <- spTransform(balleny_poly, proj4string(ice))
  ice <- crop(subset(ice, 1), extent(balleny_ice) + c(-25000, 25000, -25000, 25000))
  ice_utm <- reproject(ice, proj4string(balleny_poly_utm), program = "raster", method = "ngb")
  
  plot(ice_utm, main = paste0("0", i, "/02/2015"))
  plot(balleny_poly_utm, add = TRUE, col = "grey")
  text(SpatialPoints(coordinates(ice_utm)), round(values(ice_utm)))
  points(segdata$x[as.character(as.Date(segdata$datetime)) == paste0("2015-02-0", i)], segdata$y[as.character(as.Date(segdata$datetime)) == paste0("2015-02-0", i)], pch = 19)
  
}

#------------------------------ SEA ICE DATA AMSR2 ------------------------------#

#sea ice data from http://www.iup.uni-bremen.de:8084/amsr2data/
#converted h4 files to h5 using h4toh5convert from cmd

#get file information
gdalinfo("~/Lisa/phd/Balleny Islands/remote data/sea ice/ice.h5")
gdalinfo("~/Lisa/phd/Balleny Islands/remote data/sea ice/ice_coords.h5") #coordinates stored in separate file

#plot raster of each day
for (i in 3:6) {
  
  ice <- h5read(paste0("~/Lisa/phd/Balleny Islands/remote data/sea ice/asi-AMSR2-s6250-2015020", i, "-v5.h5"), "ASI Ice Concentration")
  ice_lats <- h5read("~/Lisa/phd/Balleny Islands/remote data/sea ice/ice_coords.h5", "Latitudes")
  ice_longs <- h5read("~/Lisa/phd/Balleny Islands/remote data/sea ice/ice_coords.h5", "Longitudes")
  
  cells <- ice_lats >= extent(balleny_poly)[3] & ice_lats <= extent(balleny_poly)[4] & 
    ice_longs >= extent(balleny_poly)[1] & ice_longs <= extent(balleny_poly)[2]
  
  ice <- ice[cells]
  ice[is.nan(ice)] <- -999 #set NaN to -999 because rasterize can't handle NA
  ice_lats <- ice_lats[cells]
  ice_longs <- ice_longs[cells]
  
  ice_sp <- SpatialPoints(coords = cbind(ice_longs, ice_lats), proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  ice_utm <- spTransform(ice_sp, CRSobj = CRS("+proj=utm +zone=58 +south +ellps=WGS84"))
  
  ice_grid <- raster(extent(gBuffer(survey.grid, width = 6500)))
  # Choose its dat_loc_utmolution (m)
  res(ice_grid) <- 6500
  
  ice_raster <- rasterize(ice_utm, ice_grid, field = ice, FUN = mean)
  
  #remove -999 values
  newVals <- values(ice_raster)
  newVals[newVals == -999] <- NA
  ice_raster <- setValues(ice_raster, newVals)
  
  plot(ice_raster, main = paste0("2015020", i), colNA = "lightgrey")
  plot(balleny_poly_utm, add = TRUE, col = "grey")
  points(segdata$x[as.character(as.Date(segdata$datetime)) == paste0("2015-02-0", i)], segdata$y[as.character(as.Date(segdata$datetime)) == paste0("2015-02-0", i)], pch = 19)
  
}




# ------------------------ PLOTS FOR PAPER ------------------------ #

plot(whale.dsm, select = 1, main = "", xlab = "Easting", ylab = "Northing", bty = "l")
plot(balleny_poly_utm, add = T, col = "grey", bty = "l")
axis(1)
axis(2)

par(mfrow = c(2, 2), oma = c(1,1,0,0), mar = c(4,4,1,1))
plot(whale.dsm, select = 2, xlab = expression(Krill~density~(gm^-2)), bty = "l")
legend(1500, 20, "(a)", bty = "n", cex = 1.5)
plot(whale.dsm, select = 3, xlab = "Salinity (ppm)", ylim = c(-10, 10), bty = "l")
legend(33.85, 13, "(b)", bty = "n", cex = 1.5)
plot(whale.dsm, select = 4, xlab = "Bottom depth (m)", ylab = "s(bottom_depth,2.23)", ylim = c(-5, 5), bty = "l")
legend(1000, 6, "(c)", bty = "n", cex = 1.5)


plot(det_function_size)


# ------------------- BATHYMETRY PLOT ------------------------- #

dsn <- "~/Lisa/phd/Balleny Islands/bathymetry/ibcso_bed_contour_500m.shp"
ogrInfo(dsn)

shape <- readShapeSpatial("~/Lisa/phd/Balleny Islands/bathymetry/ibcso_bed_contour_500m")

proj4string(shape) <- CRS("+proj=stere +lat_0=-90 +lon_0=0 +lat_ts=-65 +ellps=WGS84 +datum=WGS84 +units=m")

shape_utm <- spTransform(shape, CRS(proj4string(balleny_poly)))

shape_crop <- as(extent(balleny_poly) + c(-0.2, 0.1, -0.1, 0.2), "SpatialPolygons")
proj4string(shape_crop) <- CRS(proj4string(balleny_poly))

track_line_lat <- spTransform(track_line, CRS(proj4string(balleny_poly)))

## Clip the map
out <- gIntersection(shape_utm, shape_crop, byid=TRUE)

#plot cruise track over map

lines_ggplot <- fortify(SpatialLinesDataFrame(out[-c(31, 33, 35, 36, 41, 45, 50:55), ], data = as.data.frame(matrix(NA, nrow = 43)), match.ID = F))

ggplot(lines_ggplot, aes(x=long, y=lat, group = group)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_path(col = "grey") + 
  geom_polygon(data=fortify(balleny_poly, region="id"), aes(x=long, y=lat, group=id), color="black", fill = "grey") +
  geom_point(data = as.data.frame(coordinates(gps_lat_long)), aes(x=Longitude, y=Latitude, group = 1), col = "red") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
  annotate("text", x = 162.4, y = -66.45, label = 'atop(bold("Young"))', parse = TRUE, size = 5) +
  annotate("text", x = 163.1, y = -66.85, label = 'atop(bold("Buckle"))', parse = TRUE, size = 5) +
  annotate("text", x = 164.8, y = -67.5, label = 'atop(bold("Sturge"))', parse = TRUE, size = 5) +
  scaleBar(lon = 162, lat = -67.6, 
           distanceLon = 10, distanceLat = 5, distanceLegend = 10, 
           dist.unit = "km", orientation = FALSE)




