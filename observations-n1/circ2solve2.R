circ2solve2 <- function (p0, p1, r0, r1) {
  
  x0 <- p0[1]; y0 <- p0[2]; x1 <- p1[1]; y1 <- p1[2] # assign x and ys to individual variables
  print(p0);print(p1)
  pointB <- SpatialPoints(cbind(x0, y0))
  pointC <- SpatialPoints(cbind(x1,y1))
  distanceA2B <- r0
  distanceA2C <- r1
  
  # create the circle polygons around the points with the distances
  polyB <- gBuffer(pointB, width = distanceA2B)
  polyC <- gBuffer(pointC, width = distanceA2C)
  
  # extract the feature coordinates of the polygon - we need to intersect lines, not polygons
  coordsB <- lapply(polyB@polygons, function(x) {x@Polygons[[1]]@coords})
  coordsC <- lapply(polyC@polygons, function(x) {x@Polygons[[1]]@coords})
  
  # create circles as lines
  lineB <- SpatialLines(list(Lines(Line(coordsB[[1]]), "B")))
  lineC <- SpatialLines(list(Lines(Line(coordsC[[1]]), "C")))
  
  # find intersections
  if (is.null(gIntersection(lineB, lineC))) coords <- NULL else coords <- coordinates(gIntersection(lineB, lineC))
  print(coords)
  
  plot(lineC)
  lines(lineB)
  points(coords, col="red", pch=16)
  points(pointB, pch=4); text(coordinates(pointB) + .4, "Point B")
  points(pointC, pch=4); text(coordinates(pointC) + .4, "Point C")
  return(coords)
}

# test

p0 <- c(-1,0); p1 <- c(1,0); r0=2; r1=2
circ2solve2(p0, p1, r0, r1)
