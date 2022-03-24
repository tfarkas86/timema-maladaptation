# circ2solve is function to determine the intersecting points of
# two circles. p0 and p1 are vectors of length 2, (x,y)
# r0 and r1 are radii of the two circles

circ2solve <- function(p0, p1, r0, r1) {
  
  x0 <- p0[1]; y0 <- p0[2]; x1 <- p1[1]; y1 <- p1[2] # assign x and ys to individual variables
  dd <- sqrt(abs(y0-y1)^2 + abs(x0-x1)^2); print(dd) # distance between p0 and p1
  aa <- (r0^2 - r1^2 + dd^2) / (2*dd) # distance along p0 -> p1 to perpendicular intersection from p3
  hh <- sqrt(r1^2 - aa^2)
  p2 <- (p0 + aa * (p1 - p0)) / dd # coordinates of perpendicular intersection
  x2 <- p2[1]; y2 <- p2[2] # assign individual variables
  
  if ((dd > r0 + r1) | (dd < abs(r0-r1)) | ((dd == 0) & (r0 == r1))) {
    return("There are no solutions")
  } 
  
  if (hh == 0) {
    solution <- data.frame(x=x2, y=y2)
    row.names(solution) <- "solution"
    return(solution)
  }
  
  else {  
    
    x3.a <- (x2 + (hh * (y1 - y0))) / dd
    x3.b <- (x2 - ((hh * (y1 - y0))) / dd
    y3.a <- (y2 + ((hh * (x1 - x0))) / dd
    y3.b <- (y2 - ((hh * (x1 - x0))) / dd
    
    solutions <- data.frame(x=c(x3.a, x3.b), y=c(y3.a, y3.b))
    row.names(solutions) <- c("+ solution", "- solution")
    return(solutions) 
  }
  
}

## test

t0 <- c(-1,.5); t1 <- c(1,-.5); r0 <- 2; r1 <- 2
tt <- circ2solve(p0=t0, p1=t1, r0=r0, r1=r1)
plot(t)
