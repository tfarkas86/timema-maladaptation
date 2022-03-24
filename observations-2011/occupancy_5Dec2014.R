rm(list=ls())

dd <- read.csv("~/Dropbox/Projects/CA Arthropods/Extinction_MS/Data/obsTimData2011.csv") # load data
dd$con_mal <- ifelse(dd$host=="A", dd$conC, dd$conA)
dd$con_ad <- ifelse(dd$host=="C", dd$conC, dd$conA)
dd2 <- data.frame(rbind(dd, dd))
dd2$con <- c(dd$con_ad, dd$con_mal)
dd2$ad <- c(rep("S", nrow(dd)),rep("D", nrow(dd)))
dd$con_rat <- dd$con_mal/dd$con_ad
dd$group <- factor(rep("A", nrow(dd))) # make all observations of one group
dd$c.pcon_mal <- dd$pcon_mal-mean(dd$pcon_mal)
dd$c.con <- dd$con-mean(dd$con)
dd$c.ln.size <- dd$ln.size-mean(dd$ln.size)
dd$c.host <- ifelse(dd$host=="A", -.5, .5)
dd_long <- cnt2bin(dd, "no_mal", "no_cam")
dd$sallele <- sqrt(dd$no_str/(dd$no_str+dd$no_grn))
dd$gallele <- 1-dd$sallele
dd$pmalallele <- ifelse(dd$host=="A", dd$gallele, dd$sallele)
ad <- with(dd,data.frame(total, occupied, con_ad, con_mal, host, mal, ln.size))
write.csv(ad, "~/Dropbox/Projects/CA Arthropods/Extinction_MS/Data/obsDataShort.csv")

###### statistical analysis: drivers of maladaptation ######

# there's much more exploration to be done here, if you like. i've taken a look at some 
# obvious stuff. feel free to have a play

b_mal <- cbind(dd$no_mal, dd$no_cam) # binomal response variable for maladaptation

an1 <- lm(mal ~ con_ad + con_mal + host + ln.size + total, data=dd)
summary(an1)

an1 <- glm(b_mal ~ con_ad + con_mal + host + ln.size, family=binomial, data=dd, 
           subset=dd$total>0) # maladaptation influenced by morph frequency?
summary(an1) # actually, looks like maybe. how to check for overdispersion in this model?

comp1 <- linearHypothesis(an1, "con_ad - con_mal = 0")
comp2 <- linearHypothesis(an1, "con_ad + con_mal = 0")

an1A <- glm(b_mal ~ con_ad + con_mal + ln.size, family=binomial, data=dd, 
           subset=dd$total>0 & host=="A") # maladaptation influenced by morph frequency?
summary(an1A) # actually, looks like maybe. how to check for overdispersion in this model?

comp1A <- linearHypothesis(an1A, "con_ad - con_mal = 0")
comp2A <- linearHypothesis(an1A, "con_ad + con_mal = 0")

an1C <- glm(b_mal ~ con_ad + con_mal + ln.size, family=binomial, data=dd, 
           subset=dd$total>0 & host=="C") # maladaptation influenced by morph frequency?
summary(an1C) # actually, looks like maybe. how to check for overdispersion in this model?

comp1C <- linearHypothesis(an1C, "con_ad - con_mal = 0")
comp2C <- linearHypothesis(an1C, "con_ad + con_mal = 0")

an1D <- lm(pmalallele ~ con_ad + con_mal + ln.size, data=dd, 
            subset=dd$total>0) # maladaptation influenced by morph frequency?
summary(an1D) # actually, looks like maybe. how to check for overdispersion in this model?

comp1D <- linearHypothesis(an1D, "con_ad - con_mal = 0")
comp2D <- linearHypothesis(an1D, "con_ad + con_mal = 0")

an1.r <- glmer(b_mal ~ con_ad + con_mal + host + ln.size + (1|id), family=binomial, data=dd, subset=dd$total>0)
summary(an1.r) # random effects model accounting for resampling of bushes

comp1 <- linearHypothesis(an1.r, "con_ad - con_mal = 0")
comp2 <- linearHypothesis(an1.r, "con_ad + con_mal = 0")

an1.PQL <- glmmPQL(bin ~ 1)

an2 <- glmer(b_mal ~ c.ln.size*host + con + pcon_mal + host + total + (1|id), family=binomial, data=dd) # holding total connectivity constant
summary(an2) # doesn't seem to hold up...

an2 <- glmmPQL(bin ~ ln.size*host + con + pcon_mal + host + mal, random=~1|id, family=binomial, data=dd_long) # holding total connectivity constant
summary(an2) # doesn't seem to hold up...

an3 <- glm(b_mal ~ c.ln.size + con + pcon_mal + host + total, family=binomial, data=dd) # holding host constant constant
summary(an3) # doesn't seem to hold up...

an4 <- glmmPQL(bin ~ con_ad + con_mal + host + ln.size, random=~1|id, family=binomial, data=dd_long) # what about "maladaptive connectivity"
summary(an4) # not here either

icc <- .0001636017^2 / (.0001636017^2 + 0.9994705^2) # intra class correlation = 2.7 e-8

an5 <- glmmPQL(bin ~ pcon_mal + host + ln.size + total + con, random=~1|group, 
               correlation=corExp(form=~jitter(x)+y),
               family=binomial, data=dd_long)
summary(an5)

an6 <- glmmPQL(bin ~ pcon_mal + host, random=~1|id, 
               correlation=corExp(form=~jitter(x)+y),
               family=binomial, data=dd_long)
summary(an6)

an7 <- glm(mal ~ con_ad + con_mal + host + ln.size + total, family=gaussian, data=dd,)
summary(an7); plot(an7) #diagnostics look very bad for this model ... 

##### statistical analysis: drivers of population size ######

dd$total_1 <- dd$total-1
contrasts(dd$host) <- c(0,1)
contrasts(dd$host) <- c(1,0)
an6 <- glm(total_1 ~ con_ad + con_mal + mal + host + ln.size, data=dd[-2,], family=quasipoisson, subset=total>0)
summary(an6) 

an6 <- glm(total ~ con_ad + con_mal + host + ln.size, family=poisson, data=dd, subset=total>0)
an7 <- glm(total ~ con_ad + con_mal + host + ln.size, family=quasipoisson, data=dd, subset=total>0)
an8 <- glm(total ~ con_ad + con_mal + mal + host + ln.size, family=quasipoisson, data=dd, subset=total>0)
summary(an6)
summary(an7)
summary(an8)

err_ab <- resid(an7)
var_plot <- variogram(err_ab ~ 1, data=dd[dd$total>0,], locations= ~x+y) # variogram
plot(var_plot$dist, var_plot$gamma, # plot variogram
     xlab="distance (metres)", ylab="semivariance", main="abundance",
     ylim =c(0,2.8))
smthr <- lowess(var_plot$dist, var_plot$gamma) # make lowess smoother
lines(x=smthr$x, y=smthr$y, lwd=1, lty=1) # plot smoother


an7 <- glm.nb(total_1 ~ con_ad + con_mal + mal + host + ln.size, data=dd, subset=total>0)

an8 <- hurdle(total ~ con_ad + con_mal + ln.size + c.host, 
              data=dd, dist="poisson", zero.dist="poisson")
an9 <- zeroinfl(total ~ con_ad + con_mal + ln.size + c.host, 
                data=dd, dist="poisson")

summary(an8)

comp1 <- linearHypothesis(an8, "zero_con_ad - zero_con_mal = 0")
comp2 <- linearHypothesis(an8, "zero_con_ad + zero_con_mal = 0")

comp3 <- linearHypothesis(an8, "count_con_ad - count_con_mal = 0")
comp4 <- linearHypothesis(an8, "count_con_ad + count_con_mal = 0")
comp5 <- linearHypothesis(an6, "con_ad - con_mal = 0")

comp1 <- glht(an8, linfct=c("con_ad - con_mal == 0", "con_ad + con_mal == 0"))
summary(comp1)
comp2 <- glht(an8, linfct="con_ad - con_mal == 0")
summary(comp2)

comp2 <- glht(an6, linfct="con_ad - con_mal == 0")
summary(comp2)

ase <- .85; cse <- .78; ax <- 1.53; cx <- .40
zee <- (ax - cx) / sqrt(ase^2 + cse^2)

an7.glmm <- glmmPQL(total ~ con_ad + con_mal + ln.size + host, random=~1|group, data=dd,
                    correlation=corExp(form=~x+y), family=poisson) # add spatial autocorrelation

histogram(~ total, data=dd, type="count", breaks=seq(-.5, 30.5, by=1), 
          main="all bushes")
histogram(~ total|host, data=dd, type="count", breaks=seq(-.5, 30.5, by=1),
          main="all bushes, split by host")
histogram(~ total|cut(con_ad, 4), data=dd, type="count", breaks=seq(-.5, 30.5, by=1), 
          main="all bushes, split by conspecific connectivity")
histogram(~ total|cut(con_mal, 4), data=dd, type="count", breaks=seq(-.5, 30.5, by=1),
          main="all bushes, split by heterospecific connectivity")
histogram(~ total|cut(ln.size, 4), data=dd, type="count", breaks=seq(-.5, 30.5, by=1),
          main="all bushes, split by plant volume")

histogram(~ total, data=dd, subset=dd$total>0, type="count", 
          breaks=seq(-.5,29.5, by=1), main="bushes with 1 or more Timema")
histogram(~ total | host, data=dd, subset=dd$total>0, type="count", 
          breaks=seq(-.5,29.5, by=1), main="one or more, host")
histogram(~ total | cut(con_ad, 4), data=dd, subset=dd$total>0, type="count", 
          breaks=seq(-.5,29.5, by=1), main="one or more, conspecific connectivity")
histogram(~ total | cut(con_mal, 4), data=dd, subset=dd$total>0, type="count", 
          breaks=seq(-.5,29.5, by=1), main="one or more, heterospecific connectivity")
  histogram(~ total | cut(ln.size, 4), data=dd, subset=dd$total>0, type="count", 
            breaks=seq(-.5,29.5, by=1), main="one or more, plant volume")

hurd_err <- resid(an8)
hvar <- variogram(hurd_err ~ 1, data=dd, locations=~x+y)
plot(hvar$dist, hvar$gamma, # plot variogram
     xlab="distance (metres)", ylab="semivariance", main="occupancy")
smthr <- lowess(hvar$dist, hvar$gamma) # make lowess smoother
lines(x=smthr$x, y=smthr$y, lwd=1, lty=1) # plot smoother

# hurdle plotting?

newdata <- data.frame(c.host=0, ln.size=mean(dd$ln.size), con_ad=mean(dd$con_ad), 
                      con_mal=seq(min(dd$con_mal), max(dd$con_mal), length.out=100))
predict(an8, newdata=newdata, se.fit=TRUE)
##### statistical analysis: drivers of patch occupancy ######

an6 <- glm(occupied ~ con_ad + con_mal + ln.size + host, family=quasibinomial, data=dd)
summary(an6) 

err_occ <- resid(an7)
var_plot <- variogram(err_occ ~ 1, data=dd, locations= ~x+y) # variogram
plot(var_plot$dist, var_plot$gamma, # plot variogram
     xlab="distance (metres)", ylab="semivariance", main="occupancy",
     ylim=c(0,1))
smthr <- lowess(var_plot$dist, var_plot$gamma) # make lowess smoother
lines(x=smthr$x, y=smthr$y, lwd=1, lty=1) # plot smoother
comp <- glht(an7, linfct = c("con_ad - con_mal == 0", "con_ad + con_mal == 0"))
summary(comp)
comp <- glht(an7, linfct = "con_ad + con_mal = 0")
# compare con_mal and con_ad
x_ad <- s_an7$coefficients[2,1]
x_mal <- s_an7$coefficients[3,1]
se_ad <- s_an7$coefficients[2,2]
se_mal <- s_an7$coefficients[3,2]

x_dif <- x_ad-x_mal
se_dif <- sqrt((se_mal)^2 + (se_ad)^2)
conf_inf <- se_dif*1.96
an8 <- glm(occupied ~ host + con_ad + ln.size, family=quasibinomial, data=dd)
summary(an8)
# GAM with space, for 

# try interaction

an9 <- glmer(occupied ~ con*ad + host + ln.size + (1|id), family=binomial, data=dd2)
summary(an9)

an9 <- glm(occupied ~ con*ad + host + ln.size, family=binomial, data=dd2)
summary(an9)
some reason said not to help...?

library(mgcv) # load GAM paackage
an7.gam <- gam(occupied ~ host + con + qimm_mal + pcon_mal + ln.size, data=dd, family=binomial)
an7.gam2 <- gam(occupied ~ host + con + qimm_mal + pcon_mal + ln.size + s(x,y), data=dd, family=binomial)
summary(an7.gam)
summary(an7.gam2) # shows even more significant effect of pcon_mal (p = 0.004)

# GLMM with space!!

dd$group <- factor(rep("A", nrow(dd))) # make all observations of one group
an7.glmm <- glmmPQL(occupied ~ con_ad + con_mal + host + ln.size, random=~1|group, data=dd,
                    correlation=corExp(form=~x+y), family=binomial) # add spatial autocorrelation
an7.glmm2 <- glmmPQL(occupied ~ host+ con + qimm_mal + pcon_mal + ln.size, random=~1|group, data=dd,
                    family=binomial) # and without SAC?
summary(an7.glmm)
comp <- glht(an7.glmm, linfct = c("con_ad - con_mal == 0", "con_ad + con_mal == 0"))
summary(an7.glmm2)

anova(an7.glmm,an7.glmm2)

## pconmal

an1 <- lm(pcon_mal ~ host, data=dd)

#### plot results

# connectivities vs. occupancy
dev.off() 



dd$host2 <- ifelse(dd$host=="A", -.5, .5) 
dd$occplot <- ifelse(dd$occupied == 1, 1, 0.3)
an7 <- glm(occupied ~ con_ad + con_mal + host2 + ln.size, family=quasibinomial, data=dd)

newdata_ad <- expand.grid(host2=0, con_mal=mean(dd$con), 
                          con_ad=seq(min(dd$con_ad),max(dd$con_ad),length.out=1000), ln.size=mean(dd$ln.size))
newdata_mal <- expand.grid(host2=0, con_ad=mean(dd$con), 
                       con_mal=seq(min(dd$con_mal),max(dd$con_mal),length.out=1000), ln.size=mean(dd$ln.size))
plotdata_ad <- predict.glm(an7, newdata_ad, se.fit=TRUE, type="response")
plotdata_ad$upper <- plotdata_ad$fit+plotdata_ad$se.fit
plotdata_ad$lower <- plotdata_ad$fit-plotdata_ad$se.fit
plotdata_mal <- predict.glm(an7, newdata_mal, se.fit=TRUE, type="response")
plotdata_mal$upper <- plotdata_mal$fit+plotdata_mal$se.fit
plotdata_mal$lower <- plotdata_mal$fit-plotdata_mal$se.fit

pg_ad_x <- c(seq(min(dd$con_ad), max(dd$con_ad), length.out=1000),
             seq(max(dd$con_ad), min(dd$con_ad), length.out=1000),
             min(dd$con_ad))
pg_ad_y <- c(plotdata_ad$upper, rev(plotdata_ad$lower), plotdata_ad$upper[1])
pg_mal_x <- c(seq(min(dd$con_mal), max(dd$con_mal), length.out=1000),
              seq(max(dd$con_mal), min(dd$con_mal), length.out=1000),
              min(dd$con_mal))
pg_mal_y <- c(plotdata_mal$upper, rev(plotdata_mal$lower), plotdata_mal$upper[1])

#png("~/Dropbox/Projects/CA Arthropods/AmNatMS/Figures/ad_vs_con_occ.png")
par(mar=c(5,5,4,2))
plot(jitter(occplot, .1) ~ con_ad, data=dd, pch=1, col="red", yaxt="n", bty="o",
     xlab="host-specific connectivity", ylab="probability patch occupied", 
     xlim=c(min(dd$con_ad), max(dd$con_ad)), cex.lab=1.3, cex.axis=1.3)
points(jitter(occplot, .1) ~ con_mal, data=dd, pch=1, col="blue")
axis(2, at = c(.3, .4, .6, .8, 1), labels=c("0.0", "0.4", "0.6", "0.8", "1.0" ), cex.axis=1.3)
library(plotrix); axis.break(2, 0.35, style="slash")

points(jitter(occupied, .1) ~ con_mal, data=dd, pch=1, col="blue")
polygon(pg_ad_x, pg_ad_y, density=175, lwd=.5, col="pink")
polygon(pg_mal_x, pg_mal_y, density=175, lwd=.5, col="lightblue")
lines(newdata_ad$con_ad, plotdata_ad$fit, lwd=2, col="red")
lines(newdata_mal$con_mal, plotdata_mal$fit, lwd=2, col="blue")
#lines(newdata_mal$con_mal, plotdata_mal$upper, lwd=2, lty=2, col="blue")
#lines(newdata_mal$con_mal, plotdata_mal$lower, lwd=2, lty=2, col="blue")
#lines(newdata_ad$con_ad, plotdata_ad$upper, lwd=2, lty=2, col="red")
#lines(newdata_ad$con_ad, plotdata_ad$lower, lwd=2, lty=2, col="red")
dev.off()

### occ vs. con w/ histograms

#histograms (au:adaptive,unoccupied, hmo: maladaptive, occupied)
hau<- hist(dd$con_ad[dd$occplot==.3], breaks=seq(0, 1.1, by=.1))
hao <- hist(dd$con_ad[dd$occplot==1], breaks=seq(0, 1.1, by=.1))
hmu <- hist(dd$con_mal[dd$occplot==.3], breaks=seq(0, 1.1, by=.1))
hmo<- hist(dd$con_ad[dd$occplot==.3], breaks=seq(0, 1.1, by=.1))
hau_ct <- hau$counts
hao_ct <- hao$counts
hmu_ct <- hmu$counts
hmo_ct <- hmo$counts
scale=.1/max(c(hau_ct+hmu_ct, hao_ct+hmo_ct)) # scale to .1
hau_mids <- hau$mids
hao_mids <- hao$mids
hmu_mids <- hmu$mids
hmo_mids <- hmo$mids
hau_df <- data.frame(counts=hau_ct, mids=hau_mids)
hao_df <- data.frame(counts=hao_ct, mids=hao_mids)
hmu_df <- data.frame(counts=hmu_ct, mids=hmu_mids)
hmo_df <- data.frame(counts=hmo_ct, mids=hmo_mids)

ef <- .1 # width of histogram bars
tr <- 1.1 # "zero" for upper historgram

#### Final Occupancy by Connectivity Plot with Histograms ####

png(file="~/Dropbox/Projects/CA Arthropods/Extinction_MS/Figures/occ_con_hist.png")
par(mar=c(5,5,4,2), lheight=.8)
plot(jitter(occplot, .1) ~ con_ad, data=dd, pch=1, col="red", yaxt="n", bty="o",
     xlab="host-specific connectivity", ylab="probability patch occupied", 
     xlim=c(-.15, max(1.1)), ylim=c(.3, 1.1), type="n", las=1, xaxt="n",
     cex.lab=1.5, cex.axis=1.3)
#points(jitter(occplot, .1) ~ con_mal, data=dd, pch=1, col="blue")
axis(2, at = c(.4, .6, .8, 1), labels=c("0.4", "0.6", "0.8", "1.0" ), 
     cex.axis=1.3, las=1)
axis(1, at = c(0, .2, .4, .6, .8, 1), 
     labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cex.axis=1.3)
#library(plotrix); axis.break(2, 0.35, style="slash")


for(i in 1:nrow(hao_df)) { # conspecific connectivity occupied
  if(hao_df$counts[i] > 0)
  polygon(x=c(hao_df$mids[i] - ef/2, hao_df$mids[i] - ef/2, hao_df$mids[i] + ef/2, 
              hao_df$mids[i] + ef/2, hao_df$mids[i] - ef/2),
          y=c(tr, tr - hao_df$counts[i]*scale, tr - hao_df$counts[i]*scale, tr, tr),
          col="red", density=300, lwd=.5, border="black")
}

for(i in 1:nrow(hmo_df)) { # maladaptive connectivity occupied
  if(hmo_df$counts[i] > 0)
  polygon(x=c(hmo_df$mids[i] - ef/2, hmo_df$mids[i] - ef/2, hmo_df$mids[i] + ef/2, 
              hmo_df$mids[i] + ef/2, hmo_df$mids[i] - ef/2),
          y=c(tr - hao_df$counts[i]*scale, 
              tr - hao_df$counts[i]*scale - hmo_df$counts[i]*scale, 
              tr - hao_df$counts[i]*scale - hmo_df$counts[i]*scale, 
              tr - hao_df$counts[i]*scale, tr - hao_df$counts[i]*scale),
          col="blue", density=75, lwd=.5, border="black")
}

for(i in 1:nrow(hau_df)) { # conspecific connectivity unoccupied
  if(hau_df$counts[i] > 0)
  polygon(x=c(hau_df$mids[i] - ef/2, hau_df$mids[i] - ef/2, hau_df$mids[i] + ef/2, 
              hau_df$mids[i] + ef/2, hau_df$mids[i] - ef/2),
          y=c(.3, .3 + hau_df$counts[i]*scale, .3 + hau_df$counts[i]*scale, .3, .3), 
          col="red", density=300, lwd=.5, border="black")
}

for(i in 1:nrow(hmu_df)) { # maladaptive connectivity unoccupied
  if(hmu_df$counts[i] > 0) 
  polygon(x=c(hmu_df$mids[i] - ef/2, hmu_df$mids[i] - ef/2, hmu_df$mids[i] + ef/2, 
              hmu_df$mids[i] + ef/2, hmu_df$mids[i] - ef/2),
          y=c(.3 + hau_df$counts[i]*scale, 
              .3 + hau_df$counts[i]*scale + hmu_df$counts[i]*scale, 
              .3 + hau_df$counts[i]*scale + hmu_df$counts[i]*scale, 
              .3 + hau_df$counts[i]*scale, .3 + hau_df$counts[i]*scale),
          col="blue", density=75, lwd=.5, border="black")
}

axis(2, at=c(seq(.3, .4, length.out=3)),
     labels=c("0", "25", "50"), pos=-.02, las=1, cex.axis=1)
text(x=-.16, y=.35, labels="unoccupied", cex=1, srt=90, lheight=.8)
axis(2, at=c(seq(1, 1.1, length.out=3)),
     labels=c("50", "25", "0"), pos=-.02, las=1, cex.axis=1)
text(x=-.16, y=1.05, labels="occupied", cex=1, srt=90)

#points(jitter(occupied, .1) ~ con_mal, data=dd, pch=1, col="blue")
polygon(pg_ad_x, pg_ad_y, density=100, lwd=.5, col="red")
polygon(pg_mal_x, pg_mal_y, density=100, lwd=.5, col="blue")
lines(newdata_ad$con_ad, plotdata_ad$fit, lwd=2, col="red")
lines(newdata_mal$con_mal, plotdata_mal$fit, lwd=2, col="blue", lty=2)

dev.off()



# connectivities vs. abundance

dd$host2 <- ifelse(dd$host=="A", -.5, .5) 

an7 <- glm(total ~ con_ad + con_mal + host2 + ln.size, family=quasipoisson, data=dd, subset=total>0)

newdata_ad <- expand.grid(host2=0, con_mal=mean(dd$con), 
                          con_ad=seq(min(dd$con_ad),max(dd$con_ad),length.out=1000), ln.size=mean(dd$ln.size))
newdata_mal <- expand.grid(host2=0, con_ad=mean(dd$con), 
                           con_mal=seq(min(dd$con_mal),max(dd$con_mal),length.out=1000), ln.size=mean(dd$ln.size))
plotdata_ad <- predict.glm(an7, newdata_ad, se.fit=TRUE, type="response")
plotdata_ad$upper <- plotdata_ad$fit+plotdata_ad$se.fit
plotdata_ad$lower <- plotdata_ad$fit-plotdata_ad$se.fit
plotdata_mal <- predict.glm(an7, newdata_mal, se.fit=TRUE, type="response")
plotdata_mal$upper <- plotdata_mal$fit+plotdata_mal$se.fit
plotdata_mal$lower <- plotdata_mal$fit-plotdata_mal$se.fit

pg_ad_x <- c(seq(min(dd$con_ad), max(dd$con_ad), length.out=1000),
             seq(max(dd$con_ad), min(dd$con_ad), length.out=1000),
             min(dd$con_ad))
pg_ad_y <- c(plotdata_ad$upper, rev(plotdata_ad$lower), plotdata_ad$upper[1])
pg_mal_x <- c(seq(min(dd$con_mal), max(dd$con_mal), length.out=1000),
              seq(max(dd$con_mal), min(dd$con_mal), length.out=1000),
              min(dd$con_mal))
pg_mal_y <- c(plotdata_mal$upper, rev(plotdata_mal$lower), plotdata_mal$upper[1])

#### Final Occupancy by Connectivity Plot with Histograms ####
png(file="~/Dropbox/Projects/CA Arthropods/Extinction_MS/Figures/abun_con_hist.png")
par(mar=c(5,5,4,2), las=1)
plot(jitter(total, .1) ~ con_ad, data=dd[dd$total>0,], pch=1, col="red", bty="o",
     xlab="host-specific connectivity", ylab=expression(paste(italic('Timema '), "abundance")), 
     xlim=c(min(c(dd$con_ad, dd$con_mal)) ,max(c(dd$con_ad, dd$con_mal))), cex.lab=1.5, cex.axis=1.3)
points(total ~ con_mal, data=dd[dd$total>0,], pch=2, col="blue")
#axis(2, at = c(.3, .4, .6, .8, 1), labels=c("0.0", "0.4", "0.6", "0.8", "1.0" ), cex.lab=1.3)
#axis.break(2, 0.35, style="slash")

#points(jitter(occupied, .1) ~ con_mal, data=dd, pch=1, col="blue")
polygon(pg_ad_x, pg_ad_y, density=100, lwd=.5, col="red")
polygon(pg_mal_x, pg_mal_y, density=100, lwd=.5, col="blue")
lines(newdata_ad$con_ad, plotdata_ad$fit, lwd=2, col="red")
lines(newdata_mal$con_mal, plotdata_mal$fit, lwd=2, col="blue", lty=2)
#lines(newdata_mal$con_mal, plotdata_mal$upper, lwd=2, lty=2, col="blue")
#lines(newdata_mal$con_mal, plotdata_mal$lower, lwd=2, lty=2, col="blue")
#lines(newdata_ad$con_ad, plotdata_ad$upper, lwd=2, lty=2, col="red")
#lines(newdata_ad$con_ad, plotdata_ad$lower, lwd=2, lty=2, col="red")
dev.off()

# con vs occupancy

dd$host2 <- ifelse(dd$host=="A", -.5, .5) 

an7 <- glm(occupied ~ host2 + con + pcon_mal + ln.size, family=binomial, data=dd)


newdata <- expand.grid(host2=0, con=seq(from=min(dd$con), to=max(dd$con), length.out=1000),
                       pcon_mal=mean(dd$pcon_mal), ln.size=mean(dd$ln.size))
plotdata <- predict.glm(an7, newdata, se.fit=TRUE, type="response")
plotdata$upper <- plotdata$fit+plotdata$se.fit
plotdata$lower <- plotdata$fit-plotdata$se.fit

plot(jitter(occupied, .1) ~ con, data=dd, pch=16, 
     xlab="total connectivity", ylab="probability patch occupied")
lines(newdata$con, plotdata$fit, lwd=2)
lines(newdata$con, plotdata$upper, lwd=2, lty=2)
lines(newdata$con, plotdata$lower, lwd=2, lty=2)

## abundance by pconmal

dd$host2 <- ifelse(dd$host=="A", -.5, .5) 

an6 <- glm(total ~ c.ln.size*host2 + host2 + c.con*host2 + c.pcon_mal*host2 + mal, 
           data=dd, family=quasipoisson, subset=total>0)

newdataA <- expand.grid(host2=-.5, 
                        c.pcon_mal=seq(from=min(dd$c.pcon_mal), to=max(dd$c.pcon_mal), length.out=1000),
                        c.con=mean(dd$c.con), c.ln.size=mean(dd$c.ln.size), mal=mean(dd$mal, na.rm=TRUE))
plotdataA <- predict.glm(an6, newdataA, se.fit=TRUE, type="response")
plotdataA$upper <- plotdataA$fit+plotdataA$se.fit
plotdataA$lower <- plotdataA$fit-plotdataA$se.fit

plot(total ~ c.pcon_mal, data=dd, pch=16,subset=dd$host2==-.5, col="blue", ylim=c(0,15),
     xlab="alternate-host connectivity", ylab="Timema abundance")
lines(newdataA$c.pcon_mal, plotdataA$fit, lwd=2, col="blue")
lines(newdataA$c.pcon_mal, plotdataA$upper, lwd=2, lty=2, col="blue")
lines(newdataA$c.pcon_mal, plotdataA$lower, lwd=2, lty=2, col="blue")

newdataC <- expand.grid(host2=.5, 
                        c.pcon_mal=seq(from=min(dd$c.pcon_mal), to=max(dd$c.pcon_mal), length.out=1000),
                        c.con=mean(dd$c.con), c.ln.size=mean(dd$c.ln.size), mal=mean(dd$mal, na.rm=TRUE))
plotdataC <- predict.glm(an6, newdataC, se.fit=TRUE, type="response")
plotdataC$upper <- plotdataC$fit+plotdataC$se.fit
plotdataC$lower <- plotdataC$fit-plotdataC$se.fit

points(total ~ c.pcon_mal, data=dd, pch=16,subset=dd$host2==.5, col="orange")
lines(newdataC$c.pcon_mal, plotdataC$fit, lwd=2, col="orange")
lines(newdataC$c.pcon_mal, plotdataC$upper, lwd=2, lty=2, col="orange")
lines(newdataC$c.pcon_mal, plotdataC$lower, lwd=2, lty=2, col="orange")

## abundance by con
dd$host2 <- ifelse(dd$host=="A", -.5, .5) 

an1 <- glm(total ~ c.ln.size*host + c.con*host + c.pcon_mal*host + mal,
           data=dd, family=quasipoisson, subset=total>0)
newdataA <- expand.grid(host="A", 
                        c.con=seq(from=min(dd$c.con), to=max(dd$c.con), length.out=1000),
                        c.pcon_mal=mean(dd$c.pcon_mal), c.ln.size=mean(dd$c.ln.size), 
                        mal=mean(dd$mal, na.rm=TRUE))
plotdataA <- predict.glm(an1, newdataA, se.fit=TRUE, type="response")
plotdataA$upper <- plotdataA$fit+plotdataA$se.fit
plotdataA$lower <- plotdataA$fit-plotdataA$se.fit

plot(total ~ c.con, data=dd, pch=16,subset=dd$host2==-.5, col="blue",
     xlab="total connectivity", ylab="Timema abundance", ylim=c(0, 10), xlim=c(-.4, .6))
lines(newdataA$c.con, plotdataA$fit, lwd=2, col="blue")
lines(newdataA$c.con, plotdataA$upper, lwd=2, lty=2, col="blue")
lines(newdataA$c.con, plotdataA$lower, lwd=2, lty=2, col="blue")

newdataC <- expand.grid(host="C", 
                        c.con=seq(from=min(dd$c.con), to=max(dd$c.con), length.out=1000),
                        c.pcon_mal=mean(dd$c.pcon_mal), c.ln.size=mean(dd$c.ln.size), 
                        mal=mean(dd$mal, na.rm=TRUE))
plotdataC <- predict.glm(an1, newdataC, se.fit=TRUE, type="response")
plotdataC$upper <- plotdataC$fit+plotdataC$se.fit
plotdataC$lower <- plotdataC$fit-plotdataC$se.fit

points(total ~ c.con, data=dd, pch=16,subset=dd$host2==.5, col="orange")
lines(newdataC$c.con, plotdataC$fit, lwd=2, col="orange")
lines(newdataC$c.con, plotdataC$upper, lwd=2, lty=2, col="orange")
lines(newdataC$c.con, plotdataC$lower, lwd=2, lty=2, col="orange")

##
symbols(dd$x, dd$y, circles=abs(dd$err_occ), col=ifelse(dd$err_occ<0, "white", "black"))

### diagnostic plots ###

# occupancy

dd$err_occ <- resid(an7.glmm) # add error to data frame
locs <- list(cbind(dd$x,dd$y))
symbols(x=dd$x, y=dd$y, circles=abs(dd$err_occ)/1.2, # look at error on a map
        bg=ifelse(dd$err_occ<0, "white", "black"), inches=FALSE,
        xlab="x", ylab="y")

var_plot <- variogram(err_occ ~ 1, data=dd, locations= ~x+y) # variogram
var_plot2 <- variogram(err_occ ~ 1, data=dd, locations= ~x+y, # "" in multiple directions
                       alpha=c(0, 45, 90, 135))
plot(var_plot$dist, var_plot$gamma, # plot variogram
     xlab="distance (metres)", ylab="semivariance", main="occupancy")
smthr <- lowess(var_plot$dist, var_plot$gamma) # make lowess smoother
lines(x=smthr$x, y=smthr$y, lwd=1, lty=1) # plot smoother

# population size

dd$err_occ <- resid(an7) # add error to data frame
locs <- list(cbind(dd$x,dd$y))
symbols(x=dd$x, y=dd$y, circles=abs(dd$err_occ)/1.2, # look at error on a map
        bg=ifelse(dd$err_occ<0, "white", "black"), inches=FALSE,
        xlab="x", ylab="y")
legend("topright",c("positive", "negative"), fill=c("black", "white"),bty="n")

var_plot <- variogram(err_occ ~ 1, data=dd, locations= ~x+y) # variogram
var_plot2 <- variogram(err_occ ~ 1, data=dd, locations= ~x+y, alpha=c(0, 45, 90, 135))
plot(var_plot$dist, var_plot$gamma, ylim=c(0, .9),
     xlab="distance (metres)", ylab="semivariance")
smthr <- lowess(var_plot$dist, var_plot$gamma)
lines(x=smthr$x, y=smthr$y, lwd=1, lty=1)

# plot the map
par(bg="grey95", xaxt="n", yaxt="n", mar=c(0,0,0,0))
symbols(x=dd$x[dd$host=="A"], 
        y=dd$y[dd$host=="A"], 
        circles=dd$ln.size[dd$host=="A"]/10,
        inches=FALSE, bg=ifelse(dd$occupied[dd$host=="A"]==1, "blue", "white"),
        fg="blue", lwd=2, ylab="", xlab="")
symbols(x=dd$x[dd$host=="C"], 
        y=dd$y[dd$host=="C"], 
        circles=dd$ln.size[dd$host=="C"]/10,
        inches=FALSE, bg=ifelse(dd$occupied[dd$host=="A"]==1, "orange", "white"),
        fg="orange",add=TRUE, lwd=2)
lines(x=c(25, 35), y=c(-3, -3), lwd=2)
text(x=30, y=-4.5, labels="10 meters")


