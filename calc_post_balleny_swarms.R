#calclate post Balleny Islands swarm data based on code from Martin Cox
#date: 14/10/2016
#author: Lisa-Marie Harrison (with code from Martin Cox)


setwd("C:/Users/43439535/Documents/Lisa/phd/Balleny Islands/remote data/krill aggregations")


k38  <- read.csv("038kHzaggregationsPostBalleny.csv", header = T)
k120 <- read.csv("120kHzaggregationsPostBalleny.csv", header = T)


###Work on processed aggregations:
aggs <- k120

#Check diner (2001) aquatic living resources beam corection:
Nb <- function (Sv, L, P, Th, theta) {
  #Sv=observed Sv
  #L=observed length
  #Th=minimum threshold
  #theta=beam width
  #P=mean depth of the echotrace.
  
  dST <- Sv - Th
  
  B <- 0.44*theta * dST^0.45#detection angle
  
  B <- pi*B/180
  Nb <- L/(2*P*tan(B/2))
  return(Nb)
  
}
aggs$Nbi <- Nb(Sv = aggs$Sv_mean, L = aggs$Uncorrected_length, P = aggs$Depth_mean, Th = -70, theta = 7)
#when the variable Nbi=NaN it means the acoustic energy in the aggregation is less than the minimum data threshold.
aggs$include <- TRUE

aggs$include[is.nan(aggs$Nbi)] <- FALSE
aggs$include[aggs$Nbi < 1.5] <- FALSE #THis is a condition from the Diner 2001 paper
aggs <- aggs[aggs$include == TRUE, ]

#isolate krill based on 2 frequency dB difference:
mindB <- 0.37
maxdB <- 12

aggs$dBdiff <- k120$Sv_mean - k38$Sv_mean
aggs <- aggs[which(aggs$dBdiff >= mindB & aggs$dBdiff <= maxdB), ]
#write.csv(aggs,'AllKrillSwarmsPostBalleny.csv',row.names=FALSE)

#scale 120 kHz Sv to volumetric density
TS <- read.csv('TS.csv')

len <- as.numeric(substr(names(TS)[-1], 2, nchar(names(TS)[-1])))
#linear TS:
sigma_bs <- 10**(as.vector(unlist(TS[2, -1]))/10)
sigmaBSDF <- data.frame(len = len*1000, sigmaBS = sigma_bs)
#get trawl data:
LF <- read.csv('tan1502_biologicals_to_43.csv')

#krill only
LF   <- LF[LF$species == 'EUS', ]
KLF  <- LF$lgth*10
tKLF <- table(KLF)
sigmaBSDF$density <- 0
sigmaBSDF$density[which(sigmaBSDF$len %in% as.numeric(names(tKLF)))] <- as.vector(tKLF)/sum(tKLF)

TSbar <- 10*log10(weighted.mean(x=sigmaBSDF$sigmaBS,w=sigmaBSDF$density))
wetmassFunc <- function (l) {
  #l in mm
  #wetmass g from jarvis et al 2010 dsrII p930
  2.236e-6*l^3.314
}
TSkg <- 10*log10(1000/weighted.mean(x = wetmassFunc(sigmaBSDF$len), w = sigmaBSDF$density)) + TSbar

aggs$volDengPerm3 <- 1000*(10**((aggs$Sv_mean - TSkg)/10))
#write.csv(aggs, 'AllKrillSwarmsPostBalleny.csv', row.names = FALSE)
