segdata_new <- NULL

for (i in seq(1, nrow(segdata) - 1, by = 3)) {
  
  if (segdata$Transect.Label[i] == segdata$Transect.Label[i + 1] & segdata$Transect.Label[i] == segdata$Transect.Label[i + 2]) {
    
    temp <- colMeans(segdata[c(i:(i + 2)), ])
    temp[5] <- sum(segdata[c(i:(i + 2)), 5])
    temp[9] <- sum(segdata[c(i:(i + 2)), 9])
    temp[14] <- mean(segdata[c(i:(i + 2)), 14])
    temp[7] <- mean(segdata[c(i:(i + 2)), 7])
    
    segdata_new <- rbind(segdata_new, temp)
    
  }
  
}

segdata_new[, 7] <- 1:nrow(segdata_new)
segdata_new[, 14] <- chron(segdata_new[, 14])
names(segdata_new) <- colnames(segdata_new)
segdata_new <- as.data.frame(segdata_new)

obsdata_new <- obsdata
obsdata_new$Sample.Label <- ceiling(match(obsdata_new$Sample.Label, segdata$Sample.Label)/3)
segment.area_new <- segdata_new[, 5] * (truncation_width - 200)*2
obsdata_new <- obsdata_new[obsdata_new$distance >= 200 & obsdata_new$distance <= truncation_width, ]

obsdata_new <- na.omit(obsdata_new)

whale.dsm <- dsm(count ~ s(x, y, bs="sw", xt=list(bnd=bnd)) + log(krill) + s(density, bs = "cs") + s(bottom_depth, bs = "cs"), family = "poisson", ddf.obj = det_function, select=T, 
                 segment.data = segdata_new, observation.data = obsdata_new, method = "REML", segment.area = segment.area_new,
                 knots = soap.knots)


summary(whale.dsm)

plot(whale.dsm, pages = 1)



plot(segdata_new$density, segdata_new$number)



segdata_new$mx0_bottom_depth <- as.numeric(is.na(segdata_new$bottom_depth))
segdata_new$idx0_bottom_depth <- 1
segdata_new$idx0_bottom_depth[is.na(segdata_new$bottom_depth)] <- 2:(sum(segdata_new$mx0_bottom_depth) + 1)
segdata_new$bottom_depth[is.na(segdata_new$bottom_depth)] <- mean(na.omit(segdata_new$bottom_depth))


seg_points <- SpatialPoints(segdata_new[, c("x", "y")], proj4string = CRS(proj4string(balleny_poly_utm)))


segdata_new <- segdata_new[-which(!is.na(over(seg_points, balleny_poly_utm))), ]
