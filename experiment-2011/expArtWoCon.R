library(lme4)
library(car)
library(nlme)
library(multcomp)
library(effects)
# load and attach data

setwd("~/Dropbox/Projects/CA Arthropods/Data for R")
expDataWoCon <- read.table("artData2011WoCon.txt", header=T)
herbDataWoCon <- read.table("herbDataWoCon.txt", header=T)
expData <- cbind(expDataWoCon, pHol=herbDataWoCon$pHol)
expData$densTim5 <- (expData$count0.5-expData$count5WoT)/expData$size
expData$noTim5 <- expData$count0.5-expData$count5WoT
expData$timDiff <- (100/expData$size)-(expData$noTim5/expData$size)
expData$initTimDens <- 100/expData$size
expData$lSize <- log(expData$size)
herbData <- herbDataWoCon
lrMat <- data.frame(NA)

for (i in 1:15) {
  
  lrMat <- data.frame(block=c(lrMat$block, rep(herbData$block[i], herbData$noLvs[i])), 
                      treatment=c(as.character(lrMat$treatment), rep(as.character(herbData$treatment[i]), herbData$noLvs[i])),
                      holes=c(lrMat$holes, c(rep(1,herbData$noHoles[i]), rep(0, herbData$noLvs[i]-herbData$noHoles[i])))
  ) 
  
}


# create contrast codes

strVMalCodes <- cbind(strVMal=c(1, 1, -2), CeoVGrn=c(-1, 1, 0))                   
strVGrnCodes <- cbind(strVGrn=c(0, -1, 1), x2=c(2, -1, -1))
strVCeoCodes <- cbind(strVCeo=c(1, 0, -1), x2=c(-1, 2, -1))

# change internal contrasts

contrasts(expData$treatment) <- strVMalCodes
contrasts(lrMat$treatment) <- strVMalCodes

# perform ANOVA

an1 <- lm(densTim5 ~ treatment , data=expData)
summary(an1)
an1a <- lme(densTim5 ~ treatment,  random=~1|block, data=expData)
an2 <- lm(dens5WoT ~ treatment, data=expData)
summary(an2)
an3 <- lm(dens5Sp ~ treatment, data=expData)
summary(an3)
an4 <- lm(pHol ~ treatment, data=expData)
summary(an4)

# Poisson regression

#dummy codes

contrasts(expData$treatment) <- cbind(green=c(0,1,0), striped=c(0, 0, 1))
contrasts(expData$treatment) <- cbind(Ceanothus=c(1, 0, 0), green=c(0,1,0))

an1 <- glm(noTim5 ~ treatment + log(size), data=expData, family=quasipoisson)
summary(an1)
an1a <- glm(noTim5 ~ treatment + log(size), data=expData, family=poisson)
summary(an1a)
predAn1 <- predict(an1)
summary(lm(predAn1 ~ treatment, data=expData))
summary(glht(an1, linfct = mcp(treatment="Tukey")))
an2 <- glm(count5WoT ~ treatment + log(size), data=expData, family=poisson)
summary(an2)
summary(glht(an2, linfct = mcp(treatment="Tukey")))
an3 <- glm(rich5 ~ treatment + log(size), data=expData, family=poisson)
summary(an3)
summary(glht(an3, linfct = mcp(treatment="Tukey")))
an4 <- glm(holes ~ treatment, family=binomial, data=lrMat) # AIC = 1377.8
summary(an4)
summary(glht(an4, linfct = mcp(treatment="Tukey")))

# mixed models?

an1 <- lmer(noTim5 ~ treatment + log(size) + (1|block), data=expData, family=poisson)
summary(an1) # resDeviance = 30.73, nullDeviance = 59.87, prd = 0.487 
summary(glht(an1, linfct = mcp(treatment="Tukey")))
an1Null <- lmer(noTim5 ~ log(size) + (1|block), data=expData, family=poisson)
summary(an1Null) # deviance = 59.87
an2 <- lmer(count5WoT ~ treatment + log(size) + (1|block), data=expData, family=poisson)
summary(an2) # resDeviance = 18.61, nullDeviance = 27.83, prd = .331
an2Null <- lmer(count5WoT ~ log(size) + (1|block), data=expData, family=poisson)
summary(an2Null)
an3 <- lmer(rich5 ~ treatment + log(size) + (1|block), data=expData, family=poisson)
summary(an3) # resDeviance = 8.99, nullDeviance = 13.62, prd = 0.344
an3Null <- lmer(rich5 ~ log(size) + (1|block), data=expData, family=poisson)
summary(an3Null)
an4 <- lmer(holes ~ treatment + (1|block), data=lrMat, family=binomial)
summary(an4) # resDeviance = 1372, nullDeviance = 1394, prd = 0.0158!!?
an4Null <- lmer(holes ~ (1|block), data=lrMat, family=binomial)
summary(an4Null)
# summary stats for mixed models

timEff <- effect("treatment", an1)
summary(timEff)
artEff <- effect("treatment", an2)
summary(artEff)
richEff <- effect("treatment", an3)
summary(richEff)

# get richness x treatment in ANCOVA

mod1 <- lm(rich5 ~ treatment + log(size), data=expData)
mod2 <- lm(rich5 ~ log(size), data=expData)
richLTmnt <- residuals(lm(rich5 ~ log(size), data=expData))
mod2 <- lm(richLTmnt ~ treatment, data=expData)
plot(expData$treatment, richLTmnt)
plot(expData$treatment, predict(an2b))
richLTmnt <- richLTmnt + 5
mod2 <- glm((richLTmnt) ~ treatment, data=expData, family=poisson)

# ANOVA for Tukey's HSD

an1 <- aov(densTim5 ~ treatment, data=expData)
summary(an1)
an2 <- aov(dens5WoT ~ treatment, data=expData)
summary(an2)
an3 <- aov(dens5Sp ~ treatment, data=expData)
summary(an3)
an4 <- aov(pHol ~ treatment, data=expData)
summary(an4)

TukeyHSD(an1)
TukeyHSD(an2)
TukeyHSD(an3)
TukeyHSD(an4)

# logistic regression for herbivory

herbData <- herbDataWoCon
lrMat <- data.frame(NA)

strVMalCodes <- cbind(strVMal=c(1, 1, -2), CeoVGrn=c(-1, 1, 0))                   
strVGrnCodes <- cbind(strVGrn=c(0, -1, 1), x2=c(2, -1, -1))
strVCeoCodes <- cbind(strVCeo=c(1, 0, -1), x2=c(-1, 2, -1))
contrasts(lrMat$treatment) <- strVMalCodes
contrasts(lrMat$treatment) <- srtVGrnCodes
contrasts(lrMat$treatment) <- strVCeoCodes


for (i in 1:15) {
  
  lrMat <- data.frame(block=c(lrMat$block, rep(herbData$block[i], herbData$noLvs[i])), 
                      treatment=c(as.character(lrMat$treatment), rep(as.character(herbData$treatment[i]), herbData$noLvs[i])),
                      holes=c(lrMat$holes, c(rep(1,herbData$noHoles[i]), rep(0, herbData$noLvs[i]-herbData$noHoles[i])))
                      ) 
  
}

strVMalCodes <- cbind(strVMal=c(1, 1, -2), CeoVGrn=c(-1, 1, 0))
contrasts(lrMat$treatment) <- strVMalCodes 
contrasts(lrMat$treatment) <- strVGrnCodes
contrasts(lrMat$treatment) <- strVCeoCodes

an1 <- glm(holes ~ treatment, family=binomial, data=lrMat) # AIC = 1377.8
summary(an1)
an2 <- lmer(holes ~ treatment + (1 | block), family=binomial, data=lrMat) # AIC = 1380
summary(an2)