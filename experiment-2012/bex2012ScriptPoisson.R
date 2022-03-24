# load data and libraries

setwd("~/Dropbox/Projects/CA Arthropods/Data for R")
bexData <- read.table("bex2012data.txt", header=T)
bexHerb <- read.table("bexHerbData.txt", header=T)
herbFull <- read.table("bexHerbDataFull.txt", header=T)

morphNum <- rep(c(-1,1), 8)
cagedNum <- rep(c(1, 1, 0, 0), 4)
int <- morphNum * cagedNum

library(nlme)
library(MASS)
library(lme4)
library(glmmML)

# specify contrasts

openDummy <- cbind(c(1,0))
birdCon <- cbind(c(1,-1))
morphCon <- cbind(c(-1, 1))
closedDummy<-(cbind(c(0, 1)))

contrasts(bexData$caged) <- openDummy
contrasts(bexData$morph) <- morphCon
contrasts(bexData$caged) <- birdCon
contrasts(bexData$caged)<- closedDummy

contrasts(bexHerb$caged) <- openDummy
contrasts(bexHerb$morph) <- morphCon
contrasts(bexHerb$caged) <- closedDummy

contrasts(herbFull$caged) <- openDummy
contrasts(herbFull$morph) <- morphCon
contrasts(herbFull$caged) <- closedDummy

# random effects models

B1 <- glm(tim8 ~ morph*caged, family=poisson, data=bexData)
B2 <- glmer(tim8 ~  morph*caged + (1|block), data=bexData, 
           family=poisson)
B2c <- glmer(tim8 ~  morph + (1|block), data=bexData, 
            family=poisson, subset=bexData$caged=="C")
B2o <- glmer(tim8 ~  morph + (1|block), data=bexData, 
             family=poisson, subset=bexData$caged=="N")
B2B <- lmer(tim8 ~  morph + caged + (1|block), data=bexData, 
           family=quasipoisson)
B2L <- glmmPQL(tim8 ~ morph*caged, random=(~1|block), family=poisson, data=bexData)

timEffC <- effect("morph", B2c)
summary(timEffC)
timEffO <- effect("morph", B2o)
summary(timEffO)
artEff <- effect("treatment", an2)
summary(artEff)
richEff <- effect("treatment", an3)
summary(richEff)

# partial r2 for effect of morph in open treatment

# timema abundance
resTot <- resid(lm(tim8 ~ cagedNum + int, data=bexData))
resMorph <- resid(lm(morphNum~ cagedNum + int, data=bexData))
anPart <- lm(resTot ~ resMorph)
summary(anPart)

# arthropod abundance
resTot <- resid(lm(art ~ cagedNum + int, data=bexData))
resMorph <- resid(lm(morphNum~ cagedNum + int, data=bexData))
anPart <- lm(resTot ~ resMorph)
summary(anPart)

# arthropod richness
resTot <- resid(lm(rich ~ cagedNum + int, data=bexData))
resMorph <- resid(lm(morphNum~ cagedNum + int, data=bexData))
anPart <- lm(resTot ~ resMorph)
summary(anPart)

# herbivory
resTot <- resid(lm(pHolDam ~ cagedNum + int, data=bexHerb))
resMorph <- resid(lm(morphNum~ cagedNum + int, data=bexHerb))
anPart <- lm(resTot ~ resMorph)
summary(anPart) 

###

morpsummary(B1) # AIC = 96.03
sum <-summary(B2) # AIC = 38.11
summary(B2B) # AIC = 39.31
summary(B2L)
B4 <- glm(rich ~ morph*caged, data=bexData, family=poisson)
B5 <- lmer(rich ~ morph*caged + (1|block), family=poisson, data=bexData)

B5c <- glmer(rich ~  morph + (1|block), data=bexData, 
             family=poisson, subset=bexData$caged=="C")
B5o <- glmer(rich ~  morph + (1|block), data=bexData, 
             family=poisson, subset=bexData$caged=="N")

B5B <- lmer(rich ~ morph + caged + (1|block), family=poisson, data=bexData)

B5Q <- lmer(rich ~ morph*caged + (1|block), family=quasipoisson, data=bexData)

summary(B4) # AIC = 54.14
summary(B5) # AIC = 18.58
summary(B5B) # AIC = 19.87

richEffC <- effect("morph", B5c)
summary(richEffC)
richEffO <- effect("morph", B5o)
summary(richEffO)

B6 <- glm(art ~ morph*caged, family=poisson, data=bexData)
B7 <- lmer(art ~ morph*caged + (1|block), data=bexData, family=poisson)

B7c <- glmer(art ~  morph + (1|block), data=bexData, 
             family=poisson, subset=bexData$caged=="C")
B7o <- glmer(art ~  morph + (1|block), data=bexData, 
             family=poisson, subset=bexData$caged=="N")

B7B <- lmer(art ~ morph + caged + (1|block), family=poisson, data=bexData)

summary(B6) # AIC = 62.61
summary(B7) # AIC = 24.30
summary(B7B) # AIC = 24.49

artEffC <- effect("morph", B7c)
summary(artEffC)
artEffO <- effect("morph", B7o)
summary(artEffO)

B8 <- glm(pHoles ~ morph*caged, data=bexHerb)
B9 <- glm(pHolDam ~ morph*caged, data=bexHerb)

B10 <- glm(holDam ~ morph*caged, data=herbFull, family=binomial)

B13 <- glm(holDam ~ morph+caged, data=herbFull, family=binomial)
B11 <- lmer(holDam ~ morph*caged + (1|block), data=herbFull, family=binomial)
B11c <- glmer(holDam ~  morph + (1|block), data=herbFull, 
             family=binomial, subset=bexData$caged=="C")
B11o <- glmer(holDam ~  morph + (1|block), data=herbFull, 
             family=poisson, subset=bexData$caged=="N")
B12 <- lmer(holDam ~ morph+caged + (1|block), data=herbFull, family=binomial)

summary(B8)
summary(B9)
summary(B10)
summary(B11)
summary(B12)
summary(B13)

damEffC <- effect("morph", B11c)
summary(damEffC)
damEffO <- effect("morph", B11o)
summary(damEffO)

f1 <- (sumHoles ~ morph*caged|morph*caged)
B13 <- zeroinfl(f1, dist="poisson", link="logit", data=herbFull)
summary(B13)

## summary statsitics for mixed models

bexTimEff <- effect("morph*caged", B2)
summary(bexTimEff)
bexArtEff <- effect("morph*caged", B7)
summary(bexArtEff)
bexRichEff <- effect("morph*caged", B5)
summary(bexRichEff)
bexHerbEff <- effect("morph*caged", B11)
summary(bexHerbEff)
######### model validation

A4 <- glmmML(rich ~ morph*caged, cluster= block, family=poisson, data=bexData)
A5 <- glmmML(tim8 ~ morph*caged, cluster = block, family=poisson, data=bexData)
A6 <- glmmML(art ~ morph*caged, cluster = block, family=poisson, data=bexData)
A7 <- glmmML(pHolDam ~ morph*caged, cluster=block, data=bexHerb, family=gaussian)
summary(A4)
summary(A5)

### herbivory