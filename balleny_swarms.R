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
#probably from the same distribution ( = 0.1027)
ks.test(island_swarms$volDengPerm3, post_swarms$volDengPerm3)
ks.test(island_swarms$Corrected_length, post_swarms$Corrected_length, alternative = "greater")

#create datetime
island_swarms$datetime <- chron(dates. = paste0(substr(island_swarms$Date_S, 1, 4), "/", substr(island_swarms$Date_S, 5, 6), "/", substr(island_swarms$Date_S, 7, 8)), times. = island_swarms$Time_S, format = c(dates = "y/m/d", times = "h:m:s"), out.format = c(dates = "d/m/y", times = "h:m:s"))
post_swarms$datetime <- chron(dates. = paste0(substr(post_swarms$Date_S, 1, 4), "/", substr(post_swarms$Date_S, 5, 6), "/", substr(post_swarms$Date_S, 7, 8)), times. = post_swarms$Time_S, format = c(dates = "y/m/d", times = "h:m:s"), out.format = c(dates = "d/m/y", times = "h:m:s"))


#------------------------------- ENCOUNTER RATE --------------------------------------#

island_density <- read.csv("C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/csv/CombinedKrillDen.csv", header = T)
post_density   <- read.csv("C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/csv/post_density.csv", header = T)

island_density$Distance_vl[island_density$Distance_vl < 0] <- NA
island_distance <- (max(na.omit(island_density$Distance_vl)) - min(na.omit(island_density$Distance_vl))) * 1.852 #km

post_distance <- (max(post_density$Distance_vl) - min(post_density$Distance_vl)) * 1.852 #km

#swarms/km
nrow(island_swarms)/island_distance
nrow(post_swarms)/post_distance




#Point process on a line network

marked_points <- data.frame("longitude" = island_swarms$Lon_M, "latitude" = island_swarms$Lat_M)

#convert to utm so units are in m
xy = data.frame(marked_points$longitude, marked_points$latitude)
colnames(coordinates(xy)) <- c("lon", "lat")
proj4string(xy) <- CRS("+proj=longlat +dat_botum=WGS84")
marked_points[, 1:2] <- coordinates(spTransform(xy, CRS("+proj=utm +zone=53 ellps=WGS84")))

marked_points$tp <- (marked_points$latitude - min(marked_points$latitude))/(max(marked_points$latitude) - min(marked_points$latitude))
owin <- owin(xrange = c(min(marked_points$longitude), max(marked_points$longitude)), yrange = c(min(marked_points$latitude), max(marked_points$latitude)))
transect_bounds <- ppp(x = marked_points$longitude, y = marked_points$latitude, window = owin)
jointed_vertices <- matrix(FALSE, nrow = nrow(marked_points), ncol = nrow(marked_points))
pairs <- cbind(seq(1, nrow(marked_points) - 1), seq(2, nrow(marked_points)))
jointed_vertices[rbind(pairs, cbind(pairs[, 2], pairs[, 1]))] <- TRUE


transect_line <- linnet(vertices = transect_bounds, m = jointed_vertices)
names(marked_points) <- c("x", "y", "tp")
point_network <- lpp(X = marked_points, L = transect_line)





