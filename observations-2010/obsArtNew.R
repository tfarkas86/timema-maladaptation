data <- read.table("obsArtDataSpace.txt", header=T)

data2 <- data[ which(data$site!="LOGA"), ]
data3 <- data2[ which(data2$site!="OCA"), ]


data2$fHost <- factor(data2$host)
data3$fHost <- factor(data3$host)

op <- par(mfrow=c(4,2), mar=c(3, 3, 3, 1))
dotchart(data2$sumTim, main = "Timema abundance", group = data2$fHost)
dotchart(data2$sum5, main = "abundance > 5", group = data2$fHost)
dotchart(data2$sum10, main = "abundance > 10", group = data2$fHost)
dotchart(data2$sumRich5, main = "richness > 5", group = data2$fHost)
dotchart(data2$sumRich10, main = "richness > 10", group = data2$fHost)
dotchart(data2$mal, main = "maladaptation", group = data2$fHost)
par(op)

data2$densTim <- data2$sumTim/data2$beats
data2$LdensTim <- log(data2$densTim)
data2$dens5 <- data2$sum5/data2$beats
data2$dens10 <- data2$sum10/data2$beats
data2$rich5 <- data2$sumRich5/data2$beats
data2$rich10 <- data2$sumRich10/data2$beats

op <- par(mfrow=c(4,2), mar=c(3, 3, 3, 1))
dotchart(data2$densTim, main = "Timema density", group = data2$fHost)
dotchart(data2$LdensTim, main = "log Timema density", group = data2$fHost)
dotchart(data2$dens5, main = "density > 5", group = data2$fHost)
dotchart(data2$dens10, main = "density > 10", group = data2$fHost)
dotchart(data2$rich5, main = "density > 5", group = data2$fHost)
dotchart(data2$rich10, main = "density > 10", group = data2$fHost)
dotchart(data2$mal, main = "maladaptation", group = data2$fHost)
par(op)

f1 <- formula(dens10 ~ mal + fHost + mal*fHost)

library(nlme)

an1 <- lm(dens5 ~ mal + fHost + mal*fHost, data = data2)
an2 <- lm(dens10 ~ mal + fHost + mal*fHost, data = data2)
an3 <- lm(rich5 ~ mal + fHost + mal*fHost, data = data2)
an4 <- lm(rich10 ~ mal + fHost + mal*fHost, data = data2)
an5 <- lm(LdensTim ~ mal + fHost, data = data2)
summary(an1)
summary(an2)
summary(an3)
summary(an4)
summary(an5)

an2 <- gls(f1, correlation = corRatio(form=~ lat + lon, nugget = T), data = data2)
summary(an2)

plot(an1)
plot(an2)
plot(an3)
plot(an4)

E <- rstandard(an2)
mydata <- data.frame(E, data2$lat, data2$lon)
coordinates(mydata) <- c("data2.lat", "data2.lon")
bubble(mydata, "E", col = c("black", "grey"), main = "residuals", xlab = "x-coordinate", ylab = "y-coordinate")

Vario1 <- variogram(E ~ 1, mydata)
plot(Vario1)

Vario2 <- variogram(E ~ 1, mydata, alpha = c(0, 45, 90, 135))
plot(Vario2)

data2$lat[7] <- data2$lat[7]-.00001
data2$lat[8] <- data2$lat[8]+.00002
data2$lat[11] <- data2$lat[11]+0.00001

an1 <- lm(dens5 ~ mal + fHost + mal*fHost, data = data3)
an2 <- lm(dens10 ~ mal + fHost + mal*fHost, data = data3)
an3 <- lm(rich5 ~ mal + fHost + mal*fHost, data = data3)
an4 <- lm(rich10 ~ mal + fHost + mal*fHost, data = data3)
an5 <- lm(LdensTim ~ mal + fHost, data = data3)
summary(an1)
summary(an2)
summary(an3)
summary(an4)
summary(an5)

