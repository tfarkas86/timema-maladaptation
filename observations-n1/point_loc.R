# function locates a point on a map given the coordinates of two other points
# and the distances between all points
# point 1 is unknown, 2 and 3 are known
# distances are between 1 & 2, 1 & 3, 2 & 3
# input coords are for 2 & 3

point_loc <- function(dist, xy) {
  
  # b.index determines known point with lowest x
  # b is point from which to measure x and y distances
  # for unknown point
  b.ind <- match(min(xy$x), xy$x) 
  c.ind <- match(max(xy$x), xy$x)
  
  # determine sides based on position of B and C
  sideA <- dist[3]
  sideB <- ifelse(b.ind==1, dist[2], dist[1])
  sideC <- ifelse(b.ind==1, dist[1], dist[2])
  
  # determine irregular angles with law of cosines
  angA <- acos((sideB^2 + sideC^2 - sideA^2)/(2*sideB*sideC))
  angB <- acos((sideA^2 + sideC^2 - sideB^2)/(2*sideA*sideC))
  angC <- acos((sideB^2 + sideA^2 - sideC^2)/(2*sideB*sideA))
  
  # when yb < yc
  if (xy$y[b.ind] < xy$y[c.ind]) {
    # determine angles for right triangle along NS and ES for known points
    x.bc <- xy$x[c.ind]-xy$x[b.ind] # x distance between b and c
    y.bc <- xy$y[c.ind]-xy$y[b.ind] # y distance between b and c
    angY.BC <- atan(y.bc/x.bc) 
    angX.BC <- atan(x.bc/y.bc)
    
    # determine y dist of A from B
    
    if (angB + angY.BC == 90) {
      y.coord <- sideC + xy$y[b.ind]
      x.coord <- xy$x[b.ind] 
    }
    if (angB + angY.BC < 90) {
      y.coord <- (sideC * sin(angB + angY.BC)) + xy$y[b.ind]
      x.coord <- (sideC * cos(angB + angY.BC)) + xy$x[b.ind]
    }
    if (angB + angY.BC > 90) {
      z <- 180-(angB + angY.BC)
      y.coord <- xy$y[1] + (sideC * sin(z)) 
      x.coord <- xy$x[1] - (sideC * cos(z))
    }
  }
  else {
    x.bc <- xy$x[c.ind]-xy$x[b.ind] # x distance between b and c
    y.bc <- xy$y[b.ind]-xy$y[c.ind] # y distance between b and c
    angY.BC <- atan(y.bc/x.bc) 
    angX.BC <- atan(x.bc/y.bc)
    
    # determine y dist of A from B
    
    if (angB - angY.BC  == 90) {
      y.coord <- sideC + xy$y[b.ind]
      x.coord <- xy$x[b.ind] 
    }
    if (angB + angY.BC < 90) {
      y.coord <- (sideC * sin(angB - angY.BC)) + xy$y[b.ind]
      x.coord <- (sideC * cos(angB - angY.BC)) + xy$x[b.ind]
    }
    if (angB + angY.BC > 90) {
      z <- angB-90-angY.BC
      y.coord <- xy$y[b.ind] + (sideC * cos(z)) 
      x.coord <- xy$x[b.ind] - (sideC * sin(z))
    }
  }
  xy <- cbind(x=x.coord, y=y.coord)
  return(xy)
}
## test!

test_xy <- data.frame(x=n1$x[1:3], y=n1$y[1:3])
one2dist <- xydist(x=test_xy$x[2:1], y=test_xy$y[2:1])
one3dist <- xydist(x=test_xy$x[c(2,3)], y=test_xy$y[c(2,3)])
two3dist <- xydist(x=test_xy$x[c(1,3)], y=test_xy$y[c(1,3)])

test_dist <- c(one2dist, one3dist, two3dist)
test_cords <- test_xy[-2,]

point_loc(dist=test_dist, xy=test_cords) #success!!

test_xy


