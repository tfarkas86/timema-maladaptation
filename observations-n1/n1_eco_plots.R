####################################################################
## R script for the creation of scatter and barplots relating ######
## Timema abundance, arthropod abundance and diversity to #########
## maladaptation, connectivity, and host-plant species in #########
## the N1 metacommunity sampled in 2013.                  #########
###################################################################

rm(list=ls())

### load data, remove strange extra row, and create some new variables

n1 <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data_16Sept14.csv")
n1 <- n1[-147,] # remove strange row
n1$host <- factor(n1$host) # remove inexplicable extra factor level
n1$c.con_art5 <- n1$con_art5-mean(n1$con_art5) # centered arthropod connectivity
n1$c.mal <- n1$mal-mean(n1$mal, na.rm=TRUE) # centered maladaptation
n1_occ <- n1[n1$no_tim_gs>0,] # new data set with only plants having green or striped timema
n1$cn_rat <- n1$cn_ratio
n1$c.host <- ifelse(n1$host=="A", -.5, .5)

# maladaptation effects on Timema

an_mal <- glm(no_tim ~  mal + host + ln.vol_cub + con_tim, 
              family=quasipoisson, data=n1, subset=no_tim>0)
summary(an_mal)

# new model with numerical host

an1.1 <- glm(no_tim ~  mal + num_host + ln.vol_cub + con_tim, 
             family=quasipoisson, data=n1, subset=no_tim>0) 

# on arthropod abundance (5 mm)

an2 <- glm(no_art_5 ~ mal + ln.vol_cub + host*c.con_art5, 
           family=quasipoisson, data=n1)
summary(an2)

# on arthropod diversity (5 mm)

an3 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + host*c.con_art5, 
           family=quasipoisson, data=n1)
summary(an3)

# partial PLOTTING maladapation, connectivity, and host-species

  # timema abundance 

png(filename="~/Dropbox/Projects/CA Arthropods/N1/Figures/Thesis/tim_abundance_plots.png", 
    width=1500, height=500)
par(mfrow=c(1,3), mar=c(8,8,4,2), mgp=c(3,2,0))

# maladaptation

tim_res <- resid(glm(no_tim ~ host + ln.vol_cub + con_tim, 
                     family=quasipoisson, data=n1_occ)) # timema abundance residuals
mal_res <- resid(glm(mal ~ host+ ln.vol_cub + con_tim, data=n1_occ)) # maladaptation residuals
an1 <- glm(tim_res ~ mal_res) # residual-residual model
      # values of residual maladaptation at which to predict abdundacne
newdata2 <- data.frame(mal_res=seq(from=min(mal_res), to=max(mal_res), 
                                   length.out=50))

prd_data <- predict(an1, newdata2 , se.fit=TRUE) # predictions
plt_data <- data.frame(fit=prd_data$fit, # plotting data frame
                       upper=prd_data$fit+prd_data$se.fit, 
                       lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, 
                       mal=newdata2$mal_res) 

plot(jitter(tim_res, 2) ~ jitter(mal_res,2), pch=21, bg="grey", 
     ylab="Timema abundance", xlab="maladaptation", 
     cex=2.5, cex.axis=2.5, cex.lab=2.5, ann=FALSE, las=1)
mtext(side=1, text="maladaptation", cex=2, line=5)
mtext(side=2, text=expression(paste(italic("Timema "), "abundance")), 
      cex=2, line=5)
lines(upper ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(lower ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(fit ~ mal, data=plt_data, lwd=2)
text(-.55, 3, "A", cex=2.5)

# mal for ESA presi

tim_res <- resid(glm(no_tim ~ host + ln.vol_cub + con_tim, 
                     family=quasipoisson, data=n1_occ)) # timema abundance residuals
mal_res <- resid(glm(mal ~ host+ ln.vol_cub + con_tim, data=n1_occ)) # maladaptation residuals
an1 <- glm(mal_res ~ tim_res) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdata2 <- data.frame(tim_res=seq(from=min(tim_res), to=max(tim_res), 
                                   length.out=50))

prd_data <- predict(an1, newdata2 , se.fit=TRUE) # predictions
plt_data <- data.frame(fit=prd_data$fit, # plotting data frame
                       upper=prd_data$fit+prd_data$se.fit, 
                       lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, 
                       mal=newdata2$tim_res) 

par(mar=c(5, 6, 1, 1))

plot(jitter(mal_res, 2) ~ jitter(tim_res,2), pch=21, bg="grey", 
     ylab="Timema abundance", xlab="maladaptation", 
     cex=1.5, cex.axis=1.5, cex.lab=1.5, ann=FALSE, las=1)
mtext(side=2, text="maladaptation", cex=2, line=4)
mtext(side=1, text=expression(paste(italic("Timema "), "abundance")), 
      cex=2, line=3.5)
lines(upper ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(lower ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(fit ~ mal, data=plt_data, lwd=2)

    # connectivity
tim_res <- resid(glm(no_tim ~ host + ln.vol_cub + mal, 
                     family=quasipoisson, data=n1_occ)) # timema abundance residuals
con_res <- resid(glm(con_tim ~ host+ ln.vol_cub + mal, data=n1_occ)) # maladaptation residuals
an1 <- glm(tim_res ~ con_res) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdata2 <- data.frame(con_res=seq(from=min(con_res), to=max(con_res), length.out=50))

prd_data <- predict(an1, newdata2 , se.fit=TRUE) # predictions
plt_data <- data.frame(fit=prd_data$fit, upper=prd_data$fit+prd_data$se.fit, lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, con=newdata2$con_res) # plotting data frame


plot(jitter(tim_res, 2) ~ jitter(con_res,2), pch=21, bg="grey", xlab="connectivity", ylab="Timema abundance",
     cex=2.5, cex.axis=2.5, cex.lab=2.5, ann=FALSE)
mtext(side=1, text="connectivity", cex=2, line=5)

lines(upper ~ con, data=plt_data, type="l", lty=2, lwd=2)
lines(lower ~ con, data=plt_data, type="l", lty=2, lwd=2)
lines(fit ~ con, data=plt_data, lwd=2)
text(-2.7, 2.8, "B", cex=2.5)

    # host-plant species
newdata2 <- data.frame(host=c("A", "C"), mal=mean(n1_occ$mal), ln.vol_cub=mean(n1_occ$ln.vol_cub), 
                       con_tim=mean(n1_occ$con_tim))
prd_data <- predict(an_mal, newdata2, se.fit=TRUE)
plt_data <- data.frame(fit=prd_data$fit, upper=prd_data$fit+prd_data$se.fit, lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, host=newdata2$host)
fit <- cbind(Adenostoma=plt_data$fit[1], Ceanothus=plt_data$fit[2])
mids <- barplot(fit, ylim=c(0, max(plt_data$upper)+.1), col=c("white", "white"),
                cex=1.5, cex.axis=2.5, cex.lab=2.5, cex.names=2.5, names.arg=c(expression(italic("Adenostoma")), 
                                                                                expression(italic("Ceanothus"))))
mtext(side=1, "host-plant species", cex=2, line=5.5)
arrows(mids, plt_data$upper, mids, plt_data$lower, angle=90, code=3)
text(mids[1]/3.5, 1.45, "C", cex=2.5)

dev.off()

  # arthropod abundance 

png(filename="~/Dropbox/Projects/CA Arthropods/N1/Figures/Thesis/art_abundance_plots.png", width=1500, height=500)
par(mfrow=c(1,3), mar=c(8,8,4,2), mgp=c(3,2,0))

    # maladaptation

tim_res <- resid(glm(no_art_5 ~ host * con_art5 + ln.vol_cub, 
                     family=quasipoisson, data=n1_occ)) # timema abundance residuals
mal_res <- resid(glm(mal ~ host*con_art5 + ln.vol_cub, data=n1_occ)) # maladaptation residuals
an1 <- glm(tim_res ~ mal_res) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdata2 <- data.frame(mal_res=seq(from=min(mal_res), to=max(mal_res), length.out=50))

prd_data <- predict(an1, newdata2 , se.fit=TRUE) # predictions
plt_data <- data.frame(fit=prd_data$fit, upper=prd_data$fit+prd_data$se.fit, lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, mal=newdata2$mal_res) # plotting data frame


plot(jitter(tim_res, 2) ~ jitter(mal_res,2), pch=21, bg="grey", 
     ylab="arthropod abundance", xlab="", cex=2.5, cex.axis=2.5, cex.lab=2.5, ann=FALSE)
mtext(side=1, "maladaptation", cex=2, line=5)
mtext(side=2, "arthropod abundance", cex=2, line=5)
lines(upper ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(lower ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(fit ~ mal, data=plt_data, lwd=2)
text(-.57, 7.5, "A", cex=2.5)

#### BES Grant 2017: Art Abundance by Maladaptation ####

tim_res <- resid(glm(no_art_5 ~ host * con_art5 + ln.vol_cub, 
                     family=quasipoisson, data=n1_occ)) # timema abundance residuals
mal_res <- resid(glm(mal ~ host*con_art5 + ln.vol_cub, data=n1_occ)) # maladaptation residuals
an1 <- glm(tim_res ~ mal_res) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdata2 <- data.frame(mal_res=seq(from=min(mal_res), to=max(mal_res), length.out=50))

prd_data <- predict(an1, newdata2 , se.fit=TRUE) # predictions
plt_data <- data.frame(fit=prd_data$fit, upper=prd_data$fit+prd_data$se.fit, lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, mal=newdata2$mal_res) # plotting data frame


png(filename="~/Dropbox/Projects/CA Arthropods/N1/Figures/BES Grant/art_abund_mal.png", width=500, height=500)
par(mar=c(5,6,1,1))
plot(jitter(tim_res, 2) ~ jitter(mal_res,2), pch=21, bg="grey", 
     ylab="arthropod abundance", xlab="", cex=1.5, cex.axis=1.5, 
     cex.lab=1.5, ann=FALSE, las=1)
mtext(side=1, "maladaptation (residuals)", cex=2, line=3)
mtext(side=2, "arthropod abundance (residuals)", cex=2, line=3)
lines(upper ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(lower ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(fit ~ mal, data=plt_data, lwd=2)
text(-.55, 7.4, "A", cex=2.5)
dev.off()

    # connectivity

tim_resA <- resid(glm(no_art_5 ~ ln.vol_cub + mal, 
                     family=quasipoisson, data=n1_occ, subset=n1_occ$host=="A")) # timema abundance residuals
con_resA <- resid(glm(con_tim ~ ln.vol_cub + mal, data=n1_occ, subset=n1_occ$host=="A")) # maladaptation residuals
tim_resC <- resid(glm(no_art_5 ~ ln.vol_cub + mal, 
                      family=quasipoisson, data=n1_occ, subset=n1_occ$host=="C")) # timema abundance residuals
con_resC <- resid(glm(con_tim ~ ln.vol_cub + mal, data=n1_occ, subset=n1_occ$host=="C")) # maladaptation residuals

anA <- glm(tim_resA ~ con_resA) # residual-residual model
anC <- glm(tim_resC ~ con_resC) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdataA <- data.frame(con_resA=seq(from=min(con_resA), to=max(con_resA), length.out=50))
newdataC <- data.frame(con_resC=seq(from=min(con_resC), to=max(con_resC), length.out=50))

prd_dataA <- predict(anA, newdataA , se.fit=TRUE) # predictions
prd_dataC <- predict(anC, newdataC , se.fit=TRUE) # predictions

plt_dataA <- data.frame(fit=prd_dataA$fit, upper=prd_dataA$fit+prd_dataA$se.fit, lower=prd_dataA$fit-prd_dataA$se.fit,
                       se=prd_dataA$se.fit, con=newdataA$con_resA) # plotting data frame
plt_dataC <- data.frame(fit=prd_dataC$fit, upper=prd_dataC$fit+prd_dataC$se.fit, lower=prd_dataC$fit-prd_dataC$se.fit,
                       se=prd_dataC$se.fit, con=newdataC$con_resC) # plotting data frame


plot(jitter(tim_resA, 2) ~ jitter(con_resA,2), pch=21, bg="blue", col="blue", ann=FALSE,
     ylim=c(min(rbind(tim_resA, tim_resC)), max(rbind(tim_resA, tim_resC))), xlab="", ylab="", cex=2.5, cex.axis=2.5, cex.lab=2.5)
mtext(side=1, text="connectivity", cex=2, line=5)
points(y=jitter(tim_resC, 2), x=jitter(con_resC,2), pch=21, bg="orange", col="orange", cex=2.5, cex.axis=2.5, cex.lab=2.5)
lines(upper ~ con, data=plt_dataA, type="l", lty=2, col="blue", lwd=2)
lines(lower ~ con, data=plt_dataA, type="l", lty=2, col="blue", lwd=2)
lines(fit ~ con, data=plt_dataA, col="blue", lwd=2)
lines(upper ~ con, data=plt_dataC, type="l", lty=2, col="orange", lwd=2)
lines(lower ~ con, data=plt_dataC, type="l", lty=2, col="orange", lwd=2)
lines(fit ~ con, data=plt_dataC, col="orange", lwd=2)
# legend(1,7, legend=c("Adenostoma", "Ceanothus"), bty="n", pch=c(19, 19), col=c("blue", "orange"), cex=2)
text(-2.8, 7.5, "B", cex=2.5)

    # host-plant species

newdata2 <- data.frame(host=c("A", "C"), mal=mean(n1_occ$mal), ln.vol_cub=mean(n1_occ$ln.vol_cub), 
                       c.con_art5=mean(n1_occ$c.con_art5))
prd_data <- predict(an2, newdata2, se.fit=TRUE)
plt_data <- data.frame(fit=prd_data$fit, upper=prd_data$fit+prd_data$se.fit, lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, host=newdata2$host)
fit <- cbind(Adenostoma=plt_data$fit[1], Ceanothus=plt_data$fit[2])
mids <- barplot(fit, ylim=c(0, 3), col=c("white", "white"), xlab="", ylab="", cex=1.5, cex.axis=2.5, cex.lab=2.5, cex.names=2.5,
                names.arg=c(expression(italic("Adenostoma")), 
                            expression(italic("Ceanothus"))))
mtext(side=1, text="host-plant species", cex=2, line=5.5)
arrows(mids, plt_data$upper, mids, plt_data$lower, angle=90, code=3)
text(x=mids[1]/3.75, y=2.9, "C", cex=2.5)
dev.off()

  # arthropod species richness 

an3 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + host*c.con_art5, family=quasipoisson, data=n1)

png(filename="~/Dropbox/Projects/CA Arthropods/N1/Figures/Thesis/art_richness_plots.png", width=1500, height=500)
par(mfrow=c(1,3), mar=c(8,8,4,2), mgp=c(3,2,0))

    # maladaptation

tim_resA <- resid(glm(rich_art_5 ~ ln.vol_cub + c.con_art5, 
                      family=quasipoisson, data=n1_occ, subset=n1_occ$host=="A")) # timema abundance residuals
mal_resA <- resid(glm(mal ~ ln.vol_cub + c.con_art5, data=n1_occ, subset=n1_occ$host=="A")) # maladaptation residuals
tim_resC <- resid(glm(rich_art_5 ~ ln.vol_cub + c.con_art5, 
                      family=quasipoisson, data=n1_occ, subset=n1_occ$host=="C")) # timema abundance residuals
mal_resC <- resid(glm(mal ~ ln.vol_cub + c.con_art5, data=n1_occ, subset=n1_occ$host=="C")) # maladaptation residuals
anA <- glm(tim_resA ~ mal_resA) # residual-residual model
anC <- glm(tim_resC ~ mal_resC) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdataA <- data.frame(mal_resA=seq(from=min(mal_resA), to=max(mal_resA), length.out=50))
newdataC <- data.frame(mal_resC=seq(from=min(mal_resC), to=max(mal_resC), length.out=50))

prd_dataA <- predict(anA, newdataA , se.fit=TRUE) # predictions
prd_dataC <- predict(anC, newdataC , se.fit=TRUE) # predictions

plt_dataA <- data.frame(fit=prd_dataA$fit, upper=prd_dataA$fit+prd_dataA$se.fit, lower=prd_dataA$fit-prd_dataA$se.fit,
                        se=prd_dataA$se.fit, mal=newdataA$mal_resA) # plotting data frame
plt_dataC <- data.frame(fit=prd_dataC$fit, upper=prd_dataC$fit+prd_dataC$se.fit, lower=prd_dataC$fit-prd_dataC$se.fit,
                        se=prd_dataC$se.fit, mal=newdataC$mal_resC) # plotting data frame


plot(jitter(tim_resA, 2) ~ jitter(mal_resA,2), pch=21, bg="blue", col="blue", xlab="maladaptation",
     ylab="arthropod species richness", ylim=c(min(rbind(tim_resA, tim_resC))-.05, max(rbind(tim_resA, tim_resC))),
     cex=2.5, cex.axis=2.5, cex.lab=2.5, ann=FALSE)
points(y=jitter(tim_resC, 2), x=jitter(mal_resC,2), pch=21, bg="orange", col="orange", cex=2.5, cex.axis=2.5, cex.lab=2.5)
mtext(side=1, "maladaptation", cex=2, line=5)
mtext(side=2, "arthropod species richness", cex=2, line=5)
lines(upper ~ mal, data=plt_dataA, type="l", lty=2, col="blue", lwd=2)
lines(lower ~ mal, data=plt_dataA, type="l", lty=2, col="blue", lwd=2)
lines(fit ~ mal, data=plt_dataA, col="blue", lwd=2)
lines(upper ~ mal, data=plt_dataC, type="l", lty=2, col="orange", lwd=2)
lines(lower ~ mal, data=plt_dataC, type="l", lty=2, col="orange", lwd=2)
lines(fit ~ mal, data=plt_dataC, col="orange", lwd=2)
#legend(.4,7, legend=c("Adenostoma", "Ceanothus"), bty="n", pch=c(19, 19), col=c("blue", "orange"))

#### BES Grant 2017: Art Richness by Maladaptation ####

tim_resC <- resid(glm(rich_art_5 ~ ln.vol_cub + c.con_art5, 
                      family=quasipoisson, data=n1_occ, subset=n1_occ$host=="C")) # timema abundance residuals
mal_resC <- resid(glm(mal ~ ln.vol_cub + c.con_art5, data=n1_occ, subset=n1_occ$host=="C")) # maladaptation residuals

anC <- glm(tim_resC ~ mal_resC) # residual-residual model
# values of residual maladaptation at which to predict abdundacne

newdataC <- data.frame(mal_resC=seq(from=min(mal_resC), to=max(mal_resC), length.out=50))

prd_dataC <- predict(anC, newdataC , se.fit=TRUE) # predictions


plt_dataC <- data.frame(fit=prd_dataC$fit, upper=prd_dataC$fit+prd_dataC$se.fit, lower=prd_dataC$fit-prd_dataC$se.fit,
                        se=prd_dataC$se.fit, mal=newdataC$mal_resC) # plotting data frame

png(filename="~/Dropbox/Projects/CA Arthropods/N1/Figures/BES Grant/art_rich_mal.png", width=500, height=500)
par(mar=c(5,6,1,1))
par(mar=c(5,6,1,1))
plot(jitter(tim_resC, 2) ~ jitter(mal_resC,2), pch=21, bg="grey", 
     ylab="arthropod abundance", xlab="", cex=1.5, cex.axis=1.5, 
     cex.lab=1.5, ann=FALSE, las=1)
mtext(side=1, "maladaptation (residuals)", cex=2, line=3)
mtext(side=2, "arthropod species richness (residuals)", cex=2, line=3)
lines(upper ~ mal, data=plt_dataC, type="l", lty=2, lwd=2)
lines(lower ~ mal, data=plt_dataC, type="l", lty=2, lwd=2)
lines(fit ~ mal, data=plt_dataC, lwd=2)
text(-.55, 2.32, "B", cex=2.5)
dev.off()

    # connectivity

tim_resA <- resid(glm(rich_art_5 ~ ln.vol_cub + mal, 
                      family=quasipoisson, data=n1_occ, subset=n1_occ$host=="A")) # timema abundance residuals
con_resA <- resid(glm(con_art5 ~ ln.vol_cub + mal, data=n1_occ, subset=n1_occ$host=="A")) # maladaptation residuals
tim_resC <- resid(glm(rich_art_5 ~ ln.vol_cub + mal, 
                      family=quasipoisson, data=n1_occ, subset=n1_occ$host=="C")) # timema abundance residuals
con_resC <- resid(glm(con_art5 ~ ln.vol_cub + mal, data=n1_occ, subset=n1_occ$host=="C")) # maladaptation residuals

anA <- glm(tim_resA ~ con_resA) # residual-residual model
anC <- glm(tim_resC ~ con_resC) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdataA <- data.frame(con_resA=seq(from=min(con_resA), to=max(con_resA), length.out=50))
newdataC <- data.frame(con_resC=seq(from=min(con_resC), to=max(con_resC), length.out=50))

prd_dataA <- predict(anA, newdataA , se.fit=TRUE) # predictions
prd_dataC <- predict(anC, newdataC , se.fit=TRUE) # predictions

plt_dataA <- data.frame(fit=prd_dataA$fit, upper=prd_dataA$fit+prd_dataA$se.fit, lower=prd_dataA$fit-prd_dataA$se.fit,
                        se=prd_dataA$se.fit, con=newdataA$con_resA) # plotting data frame
plt_dataC <- data.frame(fit=prd_dataC$fit, upper=prd_dataC$fit+prd_dataC$se.fit, lower=prd_dataC$fit-prd_dataC$se.fit,
                        se=prd_dataC$se.fit, con=newdataC$con_resC) # plotting data frame


plot(jitter(tim_resA, 2) ~ jitter(con_resA,2), pch=21, bg="blue", col="blue", xlab="connectivity", ylab="", ann=FALSE,
     ylim=c(min(rbind(tim_resA, tim_resC)), max(rbind(tim_resA, tim_resC))), cex=2.5, cex.axis=2.5, cex.lab=2.5)
points(y=jitter(tim_resC, 2), x=jitter(con_resC,2), pch=21, bg="orange", col="orange", cex=2.5, cex.axis=2.5, cex.lab=2.5)
mtext(side=1, text="connectivity", cex=2, line=5)
lines(upper ~ con, data=plt_dataA, type="l", lty=2, col="blue", lwd=2)
lines(lower ~ con, data=plt_dataA, type="l", lty=2, col="blue", lwd=2)
lines(fit ~ con, data=plt_dataA, col="blue", lwd=2)
lines(upper ~ con, data=plt_dataC, type="l", lty=2, col="orange", lwd=2)
lines(lower ~ con, data=plt_dataC, type="l", lty=2, col="orange", lwd=2)
lines(fit ~ con, data=plt_dataC, col="orange", lwd=2)
#legend(2.5,7, legend=c("Adenostoma", "Ceanothus"), bty="n", pch=c(19, 19), col=c("blue", "orange"))


    # host-plant species
newdata2 <- data.frame(host=c("A", "C"), c.mal=mean(n1_occ$c.mal), ln.vol_cub=mean(n1_occ$ln.vol_cub), 
                       c.con_art5=mean(n1_occ$c.con_art5))
prd_data <- predict(an3, newdata2, se.fit=TRUE)
plt_data <- data.frame(fit=prd_data$fit, upper=prd_data$fit+prd_data$se.fit, lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, host=newdata2$host)
fit <- cbind(Adenostoma=plt_data$fit[1], Ceanothus=plt_data$fit[2])
mids <- barplot(fit, ylim=c(0, 2), col=c("white", "white"), ylab="", cex.names=2.5, cex.axis=2.5, cex.lab=2.5,
                names.arg=c(expression(italic("Adenostoma")), 
                            expression(italic("Ceanothus"))))
arrows(mids, plt_data$upper, mids, plt_data$lower, angle=90, code=3)
mtext(side=1, "host-plant species", cex=2, line=5.5)

dev.off()

# CN Ratio and maladaptation
# maladaptation

an1 <-lm(cn_rat ~ c.host + mal, data=n1)

par(mar=c(5, 5, 1, 1))

cn_res <- resid(lm(cn_rat ~ host, data=n1, subset=!is.na(mal) & !is.na(cn_rat))) # CN residuals
mal_res <- resid(lm(mal ~ host, data=n1, subset=!is.na(mal) & !is.na(cn_rat))) # maladaptation residuals
an1 <- lm(cn_res ~ mal_res) # residual-residual model
# values of residual maladaptation at which to predict abdundacne
newdata2 <- data.frame(mal_res=seq(from=min(mal_res), to=max(mal_res), 
                                   length.out=50))

prd_data <- predict(an1, newdata2 , se.fit=TRUE) # predictions
plt_data <- data.frame(fit=prd_data$fit, # plotting data frame
                       upper=prd_data$fit+prd_data$se.fit, 
                       lower=prd_data$fit-prd_data$se.fit,
                       se=prd_data$se.fit, 
                       mal=newdata2$mal_res) 


png(filename="~/Dropbox/Projects/CA Arthropods/N1/Figures/BES Grant/cn_mal.png", width=500, height=500)
par(mar=c(5,6,1,1))
plot(jitter(-cn_res, 2) ~ jitter(mal_res,2), pch=21, bg="grey", 
     cex=1.5, cex.axis=1.5, 
     cex.lab=1.5, ann=FALSE, las=1)
mtext(side=1, "maladaptation (residuals)", cex=2, line=3)
mtext(side=2, "nitrogen:carbon ratio (residuals)", cex=2, line=3)
lines(-upper ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(-lower ~ mal, data=plt_data, type="l", lty=2, lwd=2)
lines(-fit ~ mal, data=plt_data, lwd=2)
text(-.51, 3.95, "C", cex=2.5)
dev.off()





