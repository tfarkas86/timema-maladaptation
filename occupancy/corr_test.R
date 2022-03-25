

response ~ var1 + var2
var1 = a + b
var2 = a/(a+b)



dd <- vector(mode="numeric")

dd3 <- vector(mode="numeric")

# 2 variables
dd <- vector(mode="numeric")
dd2 <- vector(mode="numeric")

for (i in 1:10000) {
  
ad <- rnorm(150)
mal <- rnorm(150)
sum <- ad+mal
comp <- mal/(ad+mal)
rat <- mal/ad
dd[i] <- cor(sum, rat)
dd2[i] <- cor(sum, comp)

}

hist(dd, xlab="r", breaks=30, main=NULL)
lines(x=c(mean(dd), mean(dd)), y=c(0,7000), col="red", lty=2, lwd=2)
#text(x=0.005, y=725, labels=paste("mean",round(mean(dd),6)), col="red")
hist(dd2, xlab="r", breaks=30, main=NULL)
lines(x=c(mean(dd2), mean(dd2)), y=c(0,7000), col="red", lty=2, lwd=2)
#text(x=0.005, y=725, labels=paste("mean",round(mean(dd),6)), col="red")
# numerator variables

# constant positive
for (i in 1:10000) {
  
  aa <- rnorm(1000)
  bb <- -.0000000000000001
  sum <- aa+bb
  comp <- aa/(aa+bb)
  dd2[i] <- cor(sum, comp)
  
}

hist(dd2, xlab="r", breaks=30, main=NULL)
lines(x=c(mean(dd), mean(dd)), y=c(0,3000), col="red", lty=2, lwd=2)
text(x=0.005, y=3000, labels="mean = -0.000038", col="red")

# constant negative
for (i in 1:10000) {
  
  aa <- rnorm(1000)
  bb <- -3
  sum <- aa+bb
  comp <- aa/(aa+bb)
  dd2[i] <- cor(sum, comp)
  
}

hist(dd2, xlab="r", breaks=30, main=NULL)
lines(x=c(mean(dd), mean(dd)), y=c(0,3000), col="red", lty=2, lwd=2)
text(x=0.005, y=3000, labels="mean = -0.000038", col="red")

# denominator variable

for (i in 1:10000) {
  
  aa <- 1
  bb <- rnorm(1000)
  sum <- aa+bb
  comp <- aa/(aa+bb)
  dd2[i] <- cor(sum, comp)
  
}

hist(dd2, xlab="r", breaks=30, main=NULL)
lines(x=c(mean(dd), mean(dd)), y=c(0,3000), col="red", lty=2, lwd=2)
text(x=0.005, y=3000, labels="mean = -0.000038", col="red")

aa = bb/(aa+bb)
aa/bb = 1/(aa+bb)
aa+bb = bb/aa