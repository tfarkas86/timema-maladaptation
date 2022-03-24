setwd("~/Dropbox/Projects/CA Arthropods/Data for R")
obsDataAll <- read.table("obsTimDataAllBushes.txt", header=T)
obsData1 <- read.table("obsTimDataOcc1.txt", header=T)
qData <- read.table("newQData.txt", header=T)
qDataOcc1 <- read.table("newQDataOcc1.txt", header=T)
obsData1$bush <- 1:100
obsDataAll$propC2NC <- obsDataAll$C2NC/obsDataAll$C2N
obsDataAll$propC2NA <- obsDataAll$C2NA/obsDataAll$C2N
obsDataAll$propC2Nmal <- (obsDataAll$hostAdZ*obsDataAll$propC2NA) + 
                         (obsDataAll$hostCeZ*obsDataAll$propC2NC) 
obsDataAll$q2Nmal <- ifelse(obsDataAll$host=="A", obsDataAll$q2N, 1-obsDataAll$q2N)
obsDataAll$q2NmalA <- (obsDataAll$hostAdZ*(1-obsDataAll$q2NA)) + 
  (obsDataAll$hostCeZ*obsDataAll$q2NA) 
obsDataAll$q2NmalC <- (obsDataAll$hostAdZ*(1-obsDataAll$q2NC)) + 
  (obsDataAll$hostCeZ*obsDataAll$q2NC) 

obsData1$propC2NC <- obsData1$C2NC/obsData1$C2N
obsData1$propC2NA <- obsData1$C2NA/obsData1$C2N
obsData1$propC2Nmal <- (obsData1$hostAdZ*obsData1$propC2NA) + 
  (obsData1$hostCeZ*obsData1$propC2NC) 

obsDataAll$C2Nmal <- (obsDataAll$hostAdZ*obsDataAll$C2NA) +
                      (obsDataAll$hostCeZ*obsDataAll$C2NC)

obsData1$C2Nmal <- (obsData1$hostAdZ*obsData1$C2NA) +
  (obsData1$hostCeZ*obsData1$C2NC)

obsData1$mal_tim <- as.integer(ifelse(obsData1$host=="A", obsData1$green, obsData1$striped))
obsData1$ada_tim <- as.integer(ifelse(obsData1$host=="A", obsData1$striped, obsData1$green))
obsData1$q2Nmal <- ifelse(obsData1$host=="A", 1-obsData1$q2N, obsData1$q2N)
library(gstat)
library(nlme)
library(lme4)

##### POPULATION SIZE #######

# multiple regressions with N and 2m for pop size

obsData1$C2NCCon <- obsData1$C2NC - mean(obsData1$C2NC)

an1 <- lm(logtotal ~ mal + C2NC + C2NA + q2NA + q2NC + logsize + hostCon + C2NC*hostCon + C2NA*HostCon
          + q2NA*hostCon + q2NC*hostCon, data=obsData1)
step(an1)
an1 <- lm(logtotal ~ mal + C2NCCon + logsize + hostCon + C2NCCon*hostCon, data=obsData1)
summary(an1)
an2 <- lm(logtotal ~  mal+ hostCon + C2N + logsize, data=obsData1)
summary(an2)
an3 <- lm(logtotal ~ mal + hostCon + logsize, data=obsData1)
an4 <- lm(logtotal ~ mal + propC2Nmal + logsize + hostCon + q2NA + q2NC + C2NA + C2NC, data=obsData1)
summary(an4)

an5 <- glm(total ~ logsize + host + propC2Nmal + mal + C2N, 
           family=quasipoisson, data=obsData1, subset=obsData1$total>1)
summary(an5)

an5A <- lm(logtotal ~ logsize + propC2Nmal, data=obsData1, subset=host=="A")
summary(an5A)
an5C <- lm(logtotal ~ logsize + propC2Nmal, data=obsData1, subset=host=="C")
summary(an5C)
an6 <- lm(mal ~ hostCon + logsize + propC2Nmal, data=obsData1)
summary(an6)

f1 <- formula(logtotal ~  mal+ hostCon + C2N + logsize)
f2 <- formula(greenProNRG ~ host + q2NA + q2NC)
f3 <- formula(occupied ~ q2NC + C2N + logsize)

# plot of propC2Nmal

popData_Cmal <-expand.grid(logsize=mean(obsDataAll$logsize),propC2Nmal=seq(from=0, to=1, by=.01), 
                           hostCon=0)
popDataPred_Cmal<-data.frame(popData_Cmal,pop.fit=predict(an5,popData_Cmal,type="response", se.fit=TRUE))
popDataPred_Cmal <- data.frame(popDataPred_Cmal, popUpper=popDataPred_Cmal$pop.fit.fit+popDataPred_Cmal$pop.fit.se.fit, 
                               popLower=popDataPred_Cmal$pop.fit.fit-popDataPred_Cmal$pop.fit.se.fit)
plot(pop.fit.fit~propC2Nmal, data=popDataPred_Cmal, type="l", lwd=4, xlab="Connectivity to alternate host", ylab="ln population size", 
     cex.axis=1.5, cex.lab=1.5, ylim=c(0.4, 1.5))
lines(x=popDataPred_Cmal$propC2Nmal, y=popDataPred_Cmal$popUpper, lty=2, lwd=4)
lines(x=popDataPred_Cmal$propC2Nmal, y=popDataPred_Cmal$popLower, lty=2, lwd=4)
points(x=obsData1$propC2Nmal, y=obsData1$logtotal, pch=20)


## load required packages

install.packages('gstat')
library(gstat)
library(nlme)
library(lme4)

## make a spatial bubble plot of residuals

# extract standardized residuals from linear regression
resAn2 <- rstandard(an9) 
# combine residuals and coordinates into data frame
myData <- data.frame(resAn2, obsData1$X, obsData1$Y)
# make a "coordinates" object to pass to bubble function
coordinates(myData) <- c("obsData1.X", "obsData1.Y")
# visually evaluate autocorrelation of residuals on a map
bubble(myData, "resAn2", col=c("black", "grey")) 

## make a variogram to look for spatial autocorrelation

vario <- variogram(resAn2 ~ 1, data=myData)
plot(vario)
vario2 <- variogram(resAn2 ~ 1, data=myData, alpha=c(0,45,90,135))

B1.gls <- gls(f2, data=obsData1, family=binomial)
Vario.gls <- Variogram(B1.gls, form=~ x + y, robust=TRUE, maxDist=2000, resType="pearson")

plot(Vario.gls, smooth=T)
# adding spatial correlation

B1A.gls <- gls(f2, correlation=corSpher(form=~ X + Y, nugget=TRUE), data=obsData1)                                      
B1B.gls <- gls(f2, correlation=corLin(form=~ X + Y, nugget=TRUE), data=obsData1)
B1C.gls <- gls(f2, correlation=corRatio(form=~ X + Y, nugget=TRUE), data=obsData1)
B1D.gls <- gls(f2, correlation=corGaus(form=~ X + Y, nugget=TRUE), data=obsData1)
B1e.gls <- gls(f2, correlation=corExp(form=~ X + Y, nugget=TRUE), data=obsData1)

summary(B1.gls)
summary(B1A.gls)
summary(B1B.gls)
summary(B1C.gls)
summary(B1D.gls)
summary(B1E.gls)

##### MORPH FREQUENCY #####

# multiple regressions with N and 2m for proportion green morph

an3 <- lm(greenProNRG ~  propC2NC + q2NA + q2NC + hostCon, data=obsData1)
summary(an3)

an3 <- lm(greenProNRG ~ propC2NA, subset=obsData1$host=="C", data=obsData1)
summary(an3)

an4 <- lm(greenProNRG ~ C2N + q2N + hostCon, data=obsData1)
summary(an4)
an7 <- lm(greenProNRG ~ hostCon + q2NA + q2NC + C2N, data=obsData1) # best model
summary(an7)
an8 <- lm(mal ~ propC2Nmal, data=obsData1)
summary(an8)

#binary model for maladaptation

mal_resp <- cbind(obsData1$mal_tim, obsData1$ada_tim)
an9 <- glm(mal_resp ~ propC2Nmal + host + logsize + logtotal, family=binomial, 
           data=obsData1, subset=obsData1$total>1)
summary(an9)

resid_mal <- rstandard(an9)
myData <- data.frame(resid_mal, obsData1$X, obsData1$Y)
vario2 <- variogram(resid_mal ~ 1, data=myData)
vario2 <- variogram(resid_mal ~ 1, data=myData, alpha=c(0,45,90,135))


grn_resp <- cbind(as.integer(obsData1$green), as.integer(obsData1$striped))
an10 <- glm(grn_resp ~ propC2NC*host + logsize + logtotal, family=quasibinomial, data=obsData1)


# partial coefficients for phenotype

#host

resTot <- resid(lm(greenProNRG ~ q2NA + q2NC, data=obsData1))
resHost <- resid(lm(hostCon ~ q2NA + q2NC + C2N, data=obsData1))
anPart <- lm(resTot~resHost)
summary(anPart)

#q2NA

resTot <- resid(lm(greenProNRG ~ hostCon + q2NC + propC2NC, data=obsData1))
resQ2NA <- resid(lm(q2NA ~ hostCon + q2NC + propC2NC, data=obsData1))
anPartQA <- lm(resTot~resQ2NA)
summary(anPart)




#q2NC

resTot <- resid(lm(greenProNRG ~ hostCon + q2NA + propC2NC, data=obsData1))
resQ2NC <- resid(lm(q2NC ~ hostCon + q2NA + propC2NC, data=obsData1))
anPartQC <- lm(resTot~resQ2NC)
summary(anPart)

dev.off()
par(bg="white")
plot(resTot ~ resQ2NA, pch=16, ylab="proportion green", xlab="Gimm", col="orange", cex.lab=1.5, cex.axis=1.5)
points(y= resTot, x=resQ2NC, pch=16, col="blue")
abline(anPartQA, col="orange", lwd=4)
abline(anPartQC, col="blue", lwd=4)
legend(x=-.28, y=-.35, legend=c("Gimm-A", "Gimm-C"), col=c("orange", "blue"), pch=16, bty="n", cex=1.5)

#better figure q2NA & C
an9 <- lm(greenProNRG ~ hostCon  + q2NC + q2NA + propC2NC, data=obsData1)

phenoData_q2NA <-expand.grid(q2NA=seq(from=0, to=1, by=.01), 
                             q2NC=mean(obsData1$q2NC), propC2NC=mean(obsData1$propC2NC), 
                             hostCon=0)
phenoData_q2NC <-expand.grid(q2NC=seq(from=0, to=1, by=.01), 
                             q2NA=mean(obsData1$q2NC), propC2NC=mean(obsData1$propC2NC), 
                             hostCon=0)
phenoData_pConC <- expand.grid(proC2NC=seq(from(min(obsData1$propC2NC), to=max(obsData1$propC2NC), by=.01), 
                               q2NC=mean(obsData1$q2NC), q2NA=mean(obsData1$q2NA), hostCon=0)
                               
phenoData_q2N <- rbind(phenoData_q2NC, phenoData_q2NA)
phenoData_q2N<-data.frame(phenoData_q2N, predict.lm(an9,phenoData_q2N,type="response", se.fit=TRUE))
phenoData_q2N <- phenoData_q2N[,-c(7,8)]
phenoData_q2N <- data.frame(phenoData_q2N, q2Nhost=c(rep("C", 101), rep("A", 101)))
phenoData_q2N

phenoData_q2N <- data.frame(phenoData_q2N, 
                            phenoUpper=phenoData_q2N$fit + phenoData_q2N$se.fit, 
                            phenoLower=phenoData_q2N$fit - phenoData_q2N$se.fit)

plot(fit[phenoData_q2N$q2Nhost=="A"] ~ q2NA[phenoData_q2N$q2Nhost=="A"], 
     data=phenoData_q2N, type="n", xlab="green frequency among immigrants", ylab="frequency of green morph", 
     cex.axis=1.5, cex.lab=1.3, ylim=c(0, 1.0), xlim=c(0, 1))
points(x=obsData1$q2NC, y=obsData1$greenProNRG, pch=20, col="peachpuff")
points(x=obsData1$q2NA, y=obsData1$greenProNRG, pch=20, col="lightcyan2")
lines(y=phenoData_q2N$fit[phenoData_q2N$q2Nhost=="A"], x=phenoData_q2N$q2NA[phenoData_q2N$q2Nhost=="A"], 
     lwd=4, col="blue")
lines(x=phenoData_q2N$q2NA[phenoData_q2N$q2Nhost=="A"], 
      y=phenoData_q2N$phenoUpper[phenoData_q2N$q2Nhost=="A"], 
      lty=2, lwd=2, col="blue")
lines(x=phenoData_q2N$q2NA[phenoData_q2N$q2Nhost=="A"], 
      y=phenoData_q2N$phenoLower[phenoData_q2N$q2Nhost=="A"], 
      lty=2, lwd=2, col="blue")

lines(x=phenoData_q2N$q2NC[phenoData_q2N$q2Nhost=="C"], y=phenoData_q2N$fit[phenoData_q2N$q2Nhost=="C"],
      lty=1, lwd=4, col="darkorange3")
lines(x=phenoData_q2N$q2NC[phenoData_q2N$q2Nhost=="C"], y=phenoData_q2N$phenoUpper[phenoData_q2N$q2Nhost=="C"], 
      lty=2, lwd=2, col="darkorange3")
lines(x=phenoData_q2N$q2NC[phenoData_q2N$q2Nhost=="C"], y=phenoData_q2N$phenoLower[phenoData_q2N$q2Nhost=="C"], 
      lty=2, lwd=2, col="darkorange3")
legend(x=.05, y=1.03, legend=c("from Adenostoma, p = 0.812", "from Ceanothus,    p = 0.238"), lwd=4, col=c("blue", "darkorange3"), cex=1,
       bty="n", y.intersp=.9)

# pConC plot

phenoData_pConC <- expand.grid(propC2NC=seq(from=min(obsData1$propC2NC), to=max(obsData1$propC2NC), by=.01), 
                                           q2NC=mean(obsData1$q2NC), q2NA=mean(obsData1$q2NA), hostCon=0)
phenoData_pConC <- data.frame(phenoData_pConC, predict(an9, phenoData_pConC, type="response", se.fit=TRUE))
phenoData_pConC <- data.frame(phenoData_pConC, upper=phenoData_pConC$fit + phenoData_pConC$se.fit, 
                              lower=phenoData_pConC$fit - phenoData_pConC$se.fit)

plot(fit ~ propC2NC, data=phenoData_pConC, type="l", lwd=4, ylim=c(0, 1), xlab="Connectivity to Ceanothus", 
     ylab="frequency green morph", cex.lab=1.3, cex.axis=1.2)
lines(x=phenoData_pConC$propC2NC, y=phenoData_pConC$upper, lty=2, lwd=3)
lines(x=phenoData_pConC$propC2NC, y=phenoData_pConC$lower, lty=2, lwd=3)
points(x=obsData1$propC2NC, y=obsData1$greenProNRG, pch=16, cex=.6)
#connectivity

resTot <- resid(lm(greenProNRG ~ hostCon + q2NA + q2NC, data=obsData1))
resC2N <- resid(lm(propC2NC ~ hostCon + q2NA + q2NC, data=obsData1))
anPartpConC <- lm(resTot~resC2N)
summary(anPart)

plot(resTot ~ resC2N, pch=16, ylab="proportion green", xlab="Ceanothus connectivity", cex.lab=1.5, cex.axis=1.5)
abline(anPartpConC, lwd=4)

#### OCCUPANCY #####

an8 <- glm(occupied ~ hostAdZ + logsize + q2NA + q2NC + logsize, data=obsDataAll, family=binomial)
summary(an8)
step(an8)

an9 <- glm(occupied ~ q2NC + C2N + logsize, data=obsDataAll, family=binomial)
summary(an9)

an9 <- glm(occupied ~ hostCon + C2NA + hostCon*C2NA + logsize, data=obsDataAll, family=binomial)
summary(an9)
an10 <- glm(occupied ~ hostCeZ + q2NC + hostCeZ*q2NC + logsize, data=obsDataAll, family=binomial)
summary(an10)
an10 <- glm(occupied ~ hostAdZ + q2NA + hostAdZ*q2NA + logsize, data=obsDataAll, family=binomial)
summary(an10)
an11 <- glm(occupied ~ hostCon + logsize + C2N, data=obsDataAll, family=binomial)
summary(an11)

resTot <- resid(glm(occupied ~ hostCon + logsize + C2Nmal, family=binomial, data=obsDataAll))
resC2Nmal <- resid(glm(C2Nmal ~ hostCon + logsize, family=gaussian, data=obsDataAll))
anPart <- glm(resTot~resC2Nmal, family=binomial)
summary(anPart)

#make a datframe to hold predictions
obsDataAll_cent <- data.frame(occupied=obsDataAll$occupied, hostCon=obsDataAll$hostCon, 
                              logsize=obsDataAll$logsize-mean(obsDataAll$logsize), 
                              q2NmalA=obsDataAll$q2NmalA-mean(obsDataAll$q2NmalA), 
                              q2NmalC=obsDataAll$q2NmalC-mean(obsDataAll$q2NmalC),
                              propC2Nmal=obsDataAll$propC2Nmal-mean(obsDataAll$propC2Nmal))

an9 <- glm(occupied ~ logsize + host + propC2Nmal + C2N, data=obsDataAll, family=quasibinomial)
summary(an9)

an9A <- glm(occupied ~ logsize + propC2Nmal + q2Nmal + hostCon, data=obsDataAll, family=binomial,
            subset=host=="A")
summary(an9A)
an9C <- glm(occupied ~ logsize + propC2Nmal + q2Nmal + hostCon, data=obsDataAll, family=binomial,
            subset=host=="C")
summary(an9C)

# plot of propC2Nmal

occData_Cmal <-expand.grid(logsize=mean(obsDataAll$logsize),propC2Nmal=seq(from=0, to=1, by=.01), 
                           q2Nmal=mean(obsDataAll$q2Nmal), hostCon=0)
occDataPred_Cmal<-data.frame(occData_Cmal,occ.fit=predict(an9,occData_Cmal,type="response", se.fit=TRUE))
occDataPred_Cmal <- data.frame(occDataPred_Cmal, occUpper=occDataPred_Cmal$occ.fit.fit+occDataPred_Cmal$occ.fit.se.fit, 
                               occLower=occDataPred_Cmal$occ.fit.fit-occDataPred_Cmal$occ.fit.se.fit)
plot(occ.fit.fit~propC2Nmal, data=occDataPred_Cmal, type="l", lwd=4, xlab="Connectivity to alternate host", ylab="Probability occupied", 
     cex.axis=1.5, cex.lab=1.5, ylim=c(0, 1.0))
lines(x=occDataPred_Cmal$propC2Nmal, y=occDataPred_Cmal$occUpper, lty=2, lwd=4)
lines(x=occDataPred_Cmal$propC2Nmal, y=occDataPred_Cmal$occLower, lty=2, lwd=4)
points(x=obsDataAll$propC2Nmal, y=obsDataAll$occupied, pch=20)
 
# plot of q2Nmal


occData_qMal <-expand.grid(logsize=mean(obsDataAll$logsize),q2Nmal=seq(from=0, to=1, by=.01), 
                           propC2Nmal=mean(obsDataAll$propC2Nmal), hostCon=0)
occDataPred_qMal<-data.frame(occData_qMal,occ.fit=predict(an9,occData_qMal,type="response", se.fit=TRUE))
occDataPred_qMal <- data.frame(occDataPred_qMal, occUpper=occDataPred_qMal$occ.fit.fit+occDataPred_qMal$occ.fit.se.fit, 
                               occLower=occDataPred_qMal$occ.fit.fit-occDataPred_qMal$occ.fit.se.fit)

dev.off()
par(mar=c(rep(5, 4)))
plot(occ.fit.fit~q2Nmal, data=occDataPred_qMal, type="l", lwd=4, xlab="Connectivity to maladapted phenotype", ylab="Probability occupied", 
     cex.axis=1.5, cex.lab=1.5, ylim=c(0, 1.0))
lines(x=occDataPred_qMal$q2Nmal, y=occDataPred_qMal$occUpper, lty=2, lwd=4)
lines(x=occDataPred_qMal$q2Nmal, y=occDataPred_qMal$occLower, lty=2, lwd=4)
points(x=obsDataAll$q2Nmal, y=obsDataAll$occupied, pch=20)

setwd("~/Dropbox/Presentations/EAWAG2014/Pictures for Lucerne Talk/")
png(filename="q2Nmal_occ.png", width=500, height=500)

# partial coefficients for pop size

# maladaptation

resTot <- resid(lm(logtotal ~ logsize + C2N + hostCon, data=obsData1))
resMal <- resid(lm(mal ~ logsize + C2N + hostCon, data=obsData1))
anPart <- lm(resTot ~ resMal)
summary(anPart)
predict <- fitted(lm(logtotal ~ mal + logsize + C2N + hostCon, data=obsData1))
plot(obsData1$mal, predict)

# connectivity

resTot <- resid(lm(logtotal ~ mal + logsize + hostCon, obsData1))
resC2N <- resid(lm(C2N ~ mal + logsize + hostCon, obsData1))
anPart <- lm(resTot ~ resC2N)
summary(anPart)

# patch size

resTot <- resid(lm(logtotal ~ mal + C2N + hostCon, obsData1))
resSize <- resid(lm(logsize ~ mal + C2N + hostCon, obsData1))
anPart <- lm(resTot ~ resSize)
summary(anPart)

# host plant species

resTot <- resid(lm(logtotal ~ mal + C2N + logsize, obsData1))
resHost <- resid(lm(hostCon ~ mal + C2N + logsize, obsData1))
anPart <- lm(resTot ~ resHost)
summary(anPart)

# partial coefficients for occupancy

# q2nC

resTot <- resid(lm(occupied ~ C2N + logsize, data=obsDataAll))
resQ <- resid(lm(q2NC ~ C2N + logsize, data=obsDataAll))
anPart <- lm(resTot~resQ)
summary(anPart)

# connectivity

resTot <- resid(lm(occupied ~ q2NC + logsize, data=obsDataAll))
resCon <- resid(lm(C2N ~ q2NC + logsize, data=obsDataAll))
anPart <- lm(resTot~resCon)
summary(anPart)

# patch size

resTot <- resid(lm(occupied ~ C2N + q2NC, data=obsDataAll))
resSize <- resid(lm(logsize ~ C2N + q2NC, data=obsDataAll))
anPart <- lm(resTot~resSize)
summary(anPart)




# plotting regressions for population size 

resTot <- resid(lm(logtotal ~ mal + C2A + hostCon, obsData))
resC2A <- resid(lm(C2A ~ mal + C2A + hostCon, obsData))
anPart <- lm(resTot ~ resC2A)
summary(anPart)

plot(resMal, resTot, ann=F, pch=16)
abline(anPart, lwd=2)
title(xlab="proportion maladapted (residuals)", ylab="log population size (residuals)")

scatterplot(resTot ~ resMal | hostCon, data=obsData, smooth=F, by.groups=F, 
            boxplots=F, legend.coords="topright", grid=F)
popDataOut <- data.frame(resMal, resTot, obsData$hostCon)
write.table(popDataOut, "popTimData.txt")
read.table("popTimData.txt", header=T)

# plotting regressions for frequency of green morph

resGreen <- resid(lm(propgreenNOredgrey ~ hostCon, data=obsData))
resQ <- resid(lm(q2A ~ hostCon, data=obsData))
anPart2 <- lm(resGreen ~ resQ)
summary(anPart2)
plot(resQ, resGreen)
abline(anPart2)

greenDataOut <- data.frame(resGreen, resQ, obsData$hostCon)
write.table(greenDataOut, "greenTimData.txt")

## raw data

plot(obsData1$mal, obsData1$logtotal, type="n", xlab = "proportion maladapted", ylab = "log population size")
points(obsData1$mal[obsData1$host=="A"], obsData1$logtotal[obsData1$host=="A"], col="red")
points(obsData1$mal[obsData1$host=="C"], obsData1$logtotal[obsData1$host=="C"], col="blue")
an1 <- lm(logtotal ~ mal, data = obsData1)
abline(an1)

## wrong plot (mal vs residuals)

plot(obsData1$mal, resTot, type="n", xlab = "proportion maladapted", ylab = "resid log population size")
points(obsData1$mal[obsData1$host=="A"], resTot[obsData1$host=="A"], col="red")
points(obsData1$mal[obsData1$host=="C"], resTot[obsData1$host=="C"], col="blue")

an1 <- glm(logtotal ~  mal+ hostCon + C2N + logsize, family=poisson, data=obsData1)
summary(an1)

an2 <- lm(resTot ~ obsData1$mal)
summary(an2)

an5 <- lm(resTot ~ resMal)



an3 <- lm(obsData1$logtotal ~ resMal)
summary(an3)

an4 <- lm(mal ~ hostCon + C2N + logsize, data = obsData1)
summary(an4)

## Beckerman Plot

an1A <- glm(logtotal ~  mal + C2N + logsize + hostAdZ, family=quasipoisson, data=obsData1)
summary(an1)

newXA <- expand.grid(logsize=mean(obsData1$logsize), C2N=mean(obsData1$C2N),hostAdZ=0,
                    mal=seq(0,1,0.1))

ppA <- predict.lm(an1A, newdata=newXA, interval = "confidence")
plot.theseA <- data.frame(newXA, ppA)

an1C <- glm(logtotal ~  mal + C2N + logsize + hostCeZ, family=quasipoisson, data=obsData1)
summary(an1C)

newXC <- expand.grid(logsize=mean(obsData1$logsize), C2N=mean(obsData1$C2N),hostCeZ=0,
                     mal=seq(0,1,0.1))

ppC <- predict.lm(an1C, newdata=newXC, interval = "confidence")
plot.theseC <- data.frame(newXC, ppC)

par(mfrow=c(1,2))

plot(obsData1$mal, obsData1$logtotal, pch=19, col="lightgrey", ylim=c(-1,2))
lines(fit ~ mal, data=plot.theseA)
lines(upr ~ mal, data=plot.theseA, lty=3)
lines(lwr ~ mal, data=plot.theseA, lty=3)

# model validation

an2 <- lm(logtotal ~  mal+ hostCon + C2N + logsize  , data=obsData1)
an3 <- lm(s.logtotal ~  mal+ hostCon + C2N + logsize  , data=obsData1)


plot(an3)

s.total <-sqrt(obsData1$total)
round.total <- ceiling(obsData1$total)

an1 <- lm(total ~ mal + hostCon + C2N + logsize, data = obsData1)
an2 <- lm(s.total ~ mal + hostCon + C2N + logsize, data = obsData1)
an4 <- lm(logtotal ~ mal + hostCon + C2N + logsize, data = obsData1)

par(mfrow=c(2,4))
 plot(an1); plot(an2)

an3 <- glm(round.total ~ mal + hostCon + C2N + logsize, data = obsData1, family=poisson)
summary(an3)

# pie chart
1+
slices <- c(1/25,4/25,7/25,14/25)
labs <- c("connectivity", "host species", "maladaptation", "patch size")
par(mai=c(.2,2,.2,3))
pie(slices, labs, cex =1.3)

text(.45, .35, "16%", cex=1.1)
text(0, -.4, "56%", cex=1.1)
text(-.2,.4, "28%", cex=1.1)
text(.65,.09, "4%", cex=1.1)
