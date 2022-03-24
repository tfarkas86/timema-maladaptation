### Script for adding data to N1 ##
###################################

# eu.dist uses pythagorean theorum to determine Euclidean distance between two 
# points, x and y

eu.dist <- function (x, y) {
  dist <- sqrt((x[1] - x[2])^2 + (y[1] - y[2])^2)
  return(dist)
}

x1 <- matrix(sample(10), 5)

# function to calculate radius of Earth at given latitude (theta) in degrees

georad <- function (theta, eq.rad=6378.1370, pol.rad=6356.7523) {
  t <- theta*(pi/180) # convert theta to radians
  a <- eq.rad; b <- pol.rad # reassign arguments to simple parameters
  num1 <- (a^2 * cos(t))^2 # left-hand numerator
  num2 <- (b^2 * sin(t))^2 # right-ght numerator
  den1 <- (a * cos(t))^2   # left-hand denominator
  den2 <- (b * sin(t))^2  # right-hand denominator
  rad <- sqrt((num1+num2)/(den1+den2))
  return(rad)  
}

# function to find distance between two points in DD

geodist <- function(lon1, lon2, lat1, lat2, rad=6378.137) {
  lon1 <- lon1*(pi/180); lon2 <- lon2*(pi/180) # convert to radians
  lat1 <- lat1*(pi/180); lat2 <- lat2*(pi/180)
  dlon <- lon2 - lon1; dlat <- lat2 - lat1 # create difference vectors for lat and long
  a <- (sin(dlat/2))^2 + cos(lat1) * cos(lat2) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  dist <- rad * c
  return(dist*1000)
  }

# testing functions

lat1 <- 34.51737870; lat2 <- 34.51737924
lon1 <- -119.7964608; lon2 <- -119.7964384
radius <- georad(34.5)
geodist(lat1=lat1, lat2=lat2, lon1=lon1, lon2=lon2, rad=radius)

# get xy coordinates through UTM

n1 <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data.csv") # load raw data
n1$host <- as.factor(ifelse(n1$plant_sp=="Adenostoma", "A", "C")) # make simpler host species variable
n1$mal <- ifelse(n1$host=="A", n1$no_grn_tim/(n1$no_str_tim+n1$no_grn_tim), # calculate maladaptation
                 n1$no_grn_tim/(n1$no_str_tim+n1$no_grn_tim)) 
n1$vol_cub <- n1$length * n1$width * n1$height # calculate cubic volume
n1$vol_ell <- (4/3) * pi * n1$vol_cub # calculate ellipsoid volume
n1$ln.vol_cub <- log(n1$vol_cub); n1$ln.vol_ell <- log(n1$vol_ell) # take nat log of volumes
n1$no_tim_gs <- n1$no_str_tim + n1$no_grn_tim
n1[147,1] <- 193; n1[148,1] <- 194
n1$lat[147:148] <- c(34.51746019, 34.51742674)
n1$long[147:148] <- c(-119.7961499, -119.7961979)

n1_coords <- cbind(lon=n1$long, lat=n1$lat)
library(rgdal)
n1_utm <- project(xy=n1_coords, "+proj=utm +zone=11 +ellps=GRS80")
n1_utm[97:116,] <- NA
n1 <- data.frame(n1, utm_lon=n1_utm[,1], utm_lat=n1_utm[,2])

max_lon <- max(n1$utm_lon, na.rm=TRUE); max_lat <- max(n1$utm_lat, na.rm=TRUE)
min_lon <- min(n1$utm_lon, na.rm=TRUE); min_lat <- min(n1$utm_lat, na.rm=TRUE)
n1$x <- n1$utm_lon-min_lon; n1$y <- n1$utm_lat-min_lat

#triangulate position of unknown points using known points

known_bushes <- c(94, 90, 93, 95, 97, 83, 77, 76, 193, 194)
bush_indeces <- match(known_bushes, n1$plant_id)
n1$known <- rep(FALSE, nrow(n1))
n1$known[bush_indeces] <- rep(TRUE, length(bush_indeces))


known_coords <- data.frame(id=n1$plant_id, x=n1$x, y=n1$y)[n1$known==TRUE,]
kc <- known_coords

dist_data <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_dist_mat.csv")
dd <- dist_data

id <- 108
index<-match(id, dd$unknown)

p0 <- kc[kc$id==dd$known1[index],2:3]
p1 <- kc[kc$id==dd$known2[index],2:3]
r0 <- dd$d1[index]
r1 <- dd$d2[index]


new_coords <- circ2solve2(p0=as.matrix(p0), p1=as.matrix(p1), r0=as.vector(r0), r1=as.vector(r1))
nc <- new_coords

kc <- rbind(kc, c(id, nc[1,]))

plot(kc[,2], kc[,3], type="n")
text(kc[,2], kc[,3], labels=kc[,1])

kc <- kc[order(kc$id),]


n1$x[97:116] <- kc[9:28, 2]; n1$y[97:116] <- kc[9:28,3]


### something or other
n1 <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data.csv") # load raw data
exdata$host <- as.factor(ifelse(n1$plant_sp=="Adenostoma", "A", "C")) # make simpler host species variable
exdata$mal <- ifelse(n1$host=="A", n1$no_grn_tim/(n1$no_str_tim+n1$no_grn_tim), # calculate maladaptation
                 n1$no_grn_tim/(n1$no_str_tim+n1$no_grn_tim)) 
exdata$vol_cub <- n1$length * n1$width * n1$height # calculate cubic volume
exdata$vol_ell <- (4/3) * pi * n1$vol_cub # calculate ellipsoid volume
exdata$ln.vol_cub <- log(n1$vol_cub); n1$ln.vol_ell <- log(n1$vol_ell) # take nat log of volumes
exdata$no_tim_gs <- n1$no_str_tim + n1$no_grn_tim

exdata <- data.frame(n1$plant_id)
write.csv(exdata, file="~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/extra_n1_data.csv")

n1_coords <- cbind(lon=n1$long, lat=n1$lat)
n1_utm <- project(xy=n1_coords, "+proj=utm +zone=11 +ellps=GRS80")
n1_utm[97:116,] <- NA
n1 <- data.frame(n1, n1_utm)

## volumes


