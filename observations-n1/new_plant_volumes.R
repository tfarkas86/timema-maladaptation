# script makes predictions about plant volumes based on plants with 2 samples

rm(list=ls())

np <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/older_data/new_n1_plant_data.csv") # all new plant data
np$c.vol_1 <- np$vol_1-mean(np$vol_1) # center volume variable

nd2 <- np[1:20,] # only plants with two samples
contrasts(nd2$species) <- c(-.5, .5) # center host 

an1 <- lm(vol_0 ~ vol_1 * species, data=nd2)
summary(an1)

pp <- np[21:40, c(1,2,10)] # only plants with 2014 sample
preds <- predict(an1, nd)

par(bg="white")

plot(vol_0 ~ vol_1, data=nd2, xlab="2014 volume (cubic inches)", ylab="2013 volume (cubic inches)")
abline(an1)
lines(c(0,700000), c(0, 700000), lty=2)
points(x=pp$vol_1, y=preds, pch=19)
