# load data and libraries

setwd("~/Dropbox/Projects/CA Arthropods/Data for R")
set
bexData <- read.table("bex2012data.txt", header=T)
library(nlme)

bexData[13,4]<-4


bexData
# specify contrasts

openDummy <- cbind(c(1,0))
birdCon <- cbind(c(1,-1))
morphCon <- cbind(c(-1, 1))
contrasts(bexData$caged) <- openDummy
contrasts(bexData$morph) <- morphCon
contrasts(bexData$caged) <- birdCon

# random effects models

B1 <- gls(rich ~ 1 + morph*caged, method="REML", data=bexData)
B2 <- lme(rich ~ 1 + morph*caged, data=bexData, 
          random=~1 | block, method="REML")
B3 <- lme(rich ~ 1 + morph*caged, data=bexData, 
          random=~1 + morph|block, method="REML")

summary(B1) # AIC = 52.62
summary(B2) # AIC = 53.45
summary(B3) # AIC = 56.48

B1B <- gls(rich ~ 1 + morph*caged, method="REML", data=bexData)

summary(B1B) # AIC = 55.55

B4 <- gls(tim8 ~ 1 + morph*caged, method="REML", data=bexData)
B5 <- lme(tim8 ~ 1 + morph*caged, data=bexData, 
          random=~1 | block, method="REML")

summary(B4) # AIC = 85.56 # morph effect & interaction in expected direct: 1-tail p = 0.068 and 0.087 
summary(B5) # AIC = 87.30

B4B <- gls(tim8 ~ 1 + morph+ caged, method="REML", data=bexData)

summary(B4B) # AIC = 88.82

B6 <- gls(art ~ 1 + morph*caged, method="REML", data=bexData)
B7 <- lme(art ~ 1 + morph*caged, data=bexData, 
          random=~1 | block, method="REML")

summary(B6) # AIC = 63.95 # morph and int expected: 1-tail p= 0.0112 and 0.1204
summary(B7) # AIC = 65.19

B6B <- gls(art ~ 1 + morph + caged, method="REML", data=bexData)

summary(B6B) # AIC = 64. 88

#########

plot(bexData$caged, bexData$tim10)

nOnly <- subset(bexData, caged=="N")
cOnly <- subset(bexData, caged=="C")

summary(lm(tim12 ~ morph, data=nOnly))
summary(lm(tim12 ~ morph, data=cOnly))

summary(an2)
summary(an3)
summary(an4)

# arthropod analyes

an1 <- lm(art ~ caged, data=bexData)
an2 <- lm(art ~ morph, data=bexData)
an3 <- lm(art ~ caged + morph, data=bexData)
an4 <- lm(art ~ caged + morph + caged*morph, data=bexData)

summary(an4)

plot(bexData$morph, bexData$art)

nOnly <- subset(bexData, caged=="N")
cOnly <- subset(bexData, caged=="C")

an1 <- lm(tim8 ~ morph, data=nOnly)
summary(an1)
