####################################################
####  R script for analysis of ecological data  ####
####  from N1 Timema metacommunity. No genetics ####
####################################################

rm(list=ls())

### load data and remove strange extra row

n1 <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data_16Sept14.csv")
cnd.raw <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/cn_data.csv")

n1$cn2 <- cnd.raw$cn[match(n1$plant_id, cnd.raw$id)]
n1$pC <- cnd.raw$wC[match(n1$plant_id, cnd.raw$id)]
n1$pN <- cnd.raw$wN[match(n1$plant_id, cnd.raw$id)]
##### manipulate data and add variables #####
n1 <- n1[-147,]
n1$host <- factor(n1$host)
n1$occupied <- ifelse(n1$no_tim == 0, 0, 1)
n1$pcon_mal <- ifelse(n1$host=="A", n1$con_tim_C/n1$con_tim, n1$con_tim_A/n1$con_tim)
n1$con_mal <- ifelse(n1$host=="A", n1$con_tim_C, n1$con_tim_A)
n1$con_ad <- ifelse(n1$host=="C", n1$con_tim_C, n1$con_tim_A)
n1$c.con_art5 <- n1$con_art5-mean(n1$con_art5)
n1$c.con_art <- n1$con_art-mean(n1$con_art)
n1$c.mal <- n1$mal-mean(n1$mal, na.rm=TRUE)
n1$no_mal <- ifelse(n1$host=="C", n1$no_str_tim, n1$no_grn_tim)
n1$no_adp <- ifelse(n1$host=="A", n1$no_str_tim, n1$no_grn_tim)
n1$con_grn <- n1$con_grn_A + n1$con_grn_C 
n1$con_str <- n1$con_str_A + n1$con_str_C
n1$rcon_grn <- n1$con_grn/n1$con_tim
n1$rcon_grn_A <- n1$con_grn_A/n1$con_grn
n1$rcon_grn_C <- n1$con_grn_C/n1$con_grn
n1$con_q2N_A <- n1$q2N_A*n1$con_tim
n1$con_q2N_C <- n1$q2N_C*n1$con_tim
n1$rcon_q2N_C <- n1$q2N_C/(n1$q2N_A+n1$q2N_C)
n1$rcon_vol_C <- n1$con_vol_C/n1$con_vol
n1$rcon_tim_C <- n1$con_tim_C/n1$con_tim
n1_occ <- n1[n1$no_tim_gs>0,]
n1$dens_art_5 <- n1$no_art_5 / n1$ln.vol_cub
n1$dens_tim_5 <- n1$no_tim_5 / n1$ln.vol_cub
n1$cn_rat <- n1$cn_ratio
n1$vol_cubm <- n1$vol_cub * 0.000016
n1$any_dens <- n1$no_any/n1$vol_cubm
n1$any_dens_4 <- n1$no_any_4/n1$vol_cubm
n1$mal2 <- n1$mal^2
n1$c.mal2 <- n1$c.mal^2
# maladaptation effects on Timema

n1.2 <- n1[complete.cases(n1[,c("no_tim", "mal", "mal2", "ln.vol_cub", "cn2")]), ]

an1.1 <- glm(no_tim ~  mal + host + ln.vol_cub + con_tim, family=quasipoisson, 
           data=n1.2, subset=no_tim>4)
summary(an1.1)

an1.2 <- glm(no_tim ~  mal + host + ln.vol_cub + con_tim, family=quasipoisson, 
             data=n1, subset=no_tim>1)
summary(an1.2)

an1.3 <- glm(no_tim ~  mal + host + ln.vol_cub + con_tim + cn_rat, 
             family=quasipoisson, 
             data=n1, subset=no_tim>0)
summary(an1.3)

# final model
an1.4 <- glm(no_tim ~ c.mal + c.mal2 + ln.vol_cub + con_tim + host + cn2, 
             family=quasipoisson, data=n1, subset=no_tim > 4)
summary(an1.4)

# on arthropod abundance (5 mm)

contrasts(n1$host) <- c(1, 0)

an2.1 <- glm(no_art_5 ~ host + ln.vol_cub + host*c.con_art5 + host,
           family=quasipoisson, data=n1)
summary(an2.1)

an2.2 <- glm(no_art_5 ~ c.mal + ln.vol_cub + host * c.con_art5, 
             family=quasipoisson, data=n1, subset=no_tim>1)
summary(an2.2)

an2.3 <- glm(no_art_5 ~ mal + host + ln.vol_cub + host*c.con_art5 + cn_rat, 
             family=quasipoisson, data=n1)
summary(an2.3)

an2.4 <- glm(no_art_5 ~ mal * host + ln.vol_cub, family=quasipoisson, data=n1)
summary(an2.4)

# on biomass

an2 <- glm(biomass ~ mal, family=Gamma, start=rep(0, 2),
           data=n1[n1$biomass > 0,])
summary(an2)

# on arthropod diversity (5 mm)

contrasts(n1$host) <- c(1, 0)
an3 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + host*c.con_art5, 
           family=quasipoisson, data=n1)
summary(an3)

an3 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + host*c.con_art5, 
           family=quasipoisson, data=n1, subset=no_tim_gs > 0)
summary(an3)

contrasts(n1$host) <- c(0,1)
an3.1 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + c.con_art5, 
             family=quasipoisson, data=n1, subset=no_tim>1)
summary(an3.1)

an3.2 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + host*c.con_art5 +
               cn_rat, 
           family=quasipoisson, data=n1, subset=no_tim_gs > 3)
summary(an3.2)

## connectivity on maladaptation?

mal_nomal <- cbind(n1_occ$no_mal, n1_occ$no_adp)
grn_str <- cbind(n1_occ$no_grn_tim, n1_occ$no_str_tim)

an1 <- glm(tims ~ n1_occ$con_tim, family=binomial) # no effect of connectivity
an2 <- glm(mal ~ con_tim, data=n1_occ)

an3 <- glm(mal_nomal ~ no_art_5, family=binomial, data=n1_occ, subset=no_tim>0)
summary(an3)

an4 <- lm(mal ~ no_art_5 + ln.vol_cub + host + c.con_art5, data=n1)
summary(an4)

# maladaptation on C:N
contrasts(n1$host) <- c(0, 1)
an5 <- lm(cn_rat ~ c.mal + host, data=n1, subset=no_tim>0)
summary(an5)

an6 <- lm(cn2 ~ pC, data=n1)
an7 <- lm(cn2 ~ pN, data=n1)
summary(an6)
summary(an7)

contrasts(n1$host) <- c(0, 1)
an5 <- lm(pC ~ c.mal + host, data=n1, subset=no_tim>0)
summary(an5)

contrasts(n1$host) <- c(0, 1)
an5 <- lm(pN ~ c.mal + host + no_herb , data=n1, subset=no_tim>0)
summary(an5)

# plant growth

hist(n1$avg_growth)

an.g <- lm(avg_growth ~ host + no_art_5, data=n1)
summary(an.g)

# mass-abundance 
contrasts(n1$host) <- c(0, 1)

an.mas <- lm(slopes ~ c.mal * host, data=n1)
summary(an.mas)
### many more models ####

contrasts(n1$host) <- c(1,0)

an2 <- glm(no_art_4 ~ mal + ln.vol_cub + host*con_art4, 
           family=quasipoisson, data=n1)
an3 <- glm(no_art_5 ~ mal + ln.vol_cub + host*c.con_art5, 
           family=quasipoisson, data=n1)
an4 <- glm(no_art_6 ~ mal + ln.vol_cub + host*con_art6, 
           family=quasipoisson, data=n1)
an5 <- glm(no_art_7 ~ mal + host + ln.vol_cub + con_art7, 
           family=quasipoisson, data=n1)
an6 <- glm(no_art_8 ~ mal + host + ln.vol_cub , con_art8, 
           family=quasipoisson, data=n1)
an7 <- glm(no_art_9 ~ mal + host + ln.vol_cub + con_art9, 
           family=quasipoisson, data=n1)
an8 <- glm(no_art_10 ~ mal + host + ln.vol_cub + con_art10, 
           family=quasipoisson, data=n1)

mult_art <- cbind(n1$no_art_4, n1$no_art_5, n1$no_art_6, n1$no_art_7, 
                  n1$no_art_8, n1$no_art_9, n1$no_art_10)
an2.1 <- manova(mult_art ~ mal + host + ln.vol_cub, data=n1)
summary(an2)
summary(an3)
summary(an4)
summary(an5)
summary(an6)
summary(an7)
summary(an8)

# on arthropod diversity
contrasts(n1$host) <- c(1,0)
contrasts(n1$host) <- c(0,1)

an9 <- glm(rich_art_4 ~ host*mal + host*ln.vol_cub + host*con_art4, 
           family=quasipoisson, data=n1)
an10 <- glm(rich_art_5 ~ host*c.mal + host+ ln.vol_cub + host*c.con_art5, 
            family=quasipoisson, data=n1)
an11 <- glm(rich_art_6 ~ host*mal + host*ln.vol_cub + host*con_art6, 
            family=quasipoisson, data=n1)
an12 <- glm(rich_art_7 ~ mal + host + ln.vol_cub + con_art7, 
            family=quasipoisson, data=n1)
an13 <- glm(rich_art_8 ~ mal + host + ln.vol_cub , con_art8, 
            family=quasipoisson, data=n1)
an14 <- glm(rich_art_9 ~ mal + host + ln.vol_cub + con_art9, 
            family=quasipoisson, data=n1)
an15 <- glm(rich_art_10 ~ mal + host + ln.vol_cub + con_art10, 
            family=quasipoisson, data=n1)

summary(an9)
summary(an10)
summary(an11)
summary(an12)
summary(an13)
summary(an14)
summary(an15)

# on herbviores 

an16 <- glm(no_art_herb_4 ~ mal + host + ln.vol_cub + con_art_herb4, 
            family=quasipoisson, data=n1)
an17 <- glm(no_art_herb_5 ~ mal + host + ln.vol_cub + con_art_herb5, 
            family=quasipoisson, data=n1)
an18 <- glm(no_art_herb_6 ~ mal + host + ln.vol_cub, family=quasipoisson, data=n1)
an19 <- glm(no_art_herb_7 ~ mal + host + ln.vol_cub + con_art_herb7, 
            family=quasipoisson, data=n1)
an20 <- glm(no_art_herb_8 ~ mal + host + ln.vol_cub , con_art_herb8, 
            family=quasipoisson, data=n1)
an21 <- glm(no_art_herb_9 ~ mal + host + ln.vol_cub + con_art_herb9, 
            family=quasipoisson, data=n1)
an22 <- glm(no_art_herb_10 ~ mal + host + ln.vol_cub + con_art_herb10, 
            family=quasipoisson, data=n1)

summary(an16)
summary(an17)
summary(an18)
summary(an19)
summary(an20)
summary(an21)
summary(an22)

### Occupancy models ###

an1 <- glm(occupied ~ ln.vol_cub + pcon_mal*host + host + con_tim ,
           family=quasibinomial, data=n1)
summary(an1)

an2 <- glm(no_art_4 ~ ln.vol_cub + con_art4 + host + mal + pcon_mal*host, 
           family=quasipoisson, data=n1)
summary(an2)

## occupancy plot

par(bg="lightgrey")
with(dd, symbols(x, y, circles=ln.size/10, inches=FALSE, 
                 fg=ifelse(host=="A", "blue", "orange"), 
                 bg=ifelse(occupied, ifelse(host=="A", "blue", "orange"), "white")))
polygon(x=c(20, 20, 30, 30, 20), y=c(0, .15, .15, 0, 0), 
        density=100, col="black")
text(x=25, y=-1.3, labels="10 meters", cex=.8)

##### path analysis ######
library(sem)

model1 <- 'biomass ~ ln.vol_cub + con_art5 + host
mal ~ biomass
rich_art_5 ~ biomass '

model2 <- 'biomass ~ ln.vol_cub + con_art5 + host + mal
rich_art_5 ~ biomass + host'

f1 <- sem(model1, data=alld)
summary(f1)

# piecewise SEM
library(piecewiseSEM)

model_list1 <- list(
  glm(rich_art_5 ~ no_art_5 + ln.vol_cub + host, 
      family=quasipoisson, data=n1),
  glm(no_art_5 ~ ln.vol_cub + con_art5 + host + mal, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + con_tim + host + mal, 
      family=quasipoisson, data=n1)
)
sem.fit(model_list1, n1)
sem.coefs(model_list1, n1)
#sem.lavaan(model_list, n1)

model_list2 <- list(
  glm(no_art_5 ~ ln.vol_cub + con_art5 + host, family=quasipoisson, data=n1),
  lm(mal ~ host + no_art_5, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + con_tim + host + mal, family=quasipoisson, data=n1)
)

sem.fit(model_list2, n1)
sem.coefs(model_list2, n1)
#sem.lavaan(model_list2, n1)


an1 <- glm(dens_art_5 ~ host*c.con_art5 + mal, data=n1)
an2 <- glm(no_art_5 ~ host + mal, family=qu)

model_list3 <- list(
  glm(rich_art_5 ~ host*c.con_art5 + dens_art_5, family=quasipoisson, data=n1),
  lm(dens_art_5 ~ host*c.con_art5 + mal, data=n1)
)
sem.fit(model_list3, n1)
sem.coefs(model_list3, n1)

model_list4 <- list(
  glm(rich_art_5 ~ host*c.con_art5 + dens_art_5, family=quasipoisson, data=n1),
  lm(dens_art_5 ~ host*c.con_art5, data=n1),
  lm(mal ~ dens_art_5 + host, data=n1)
)
sem.fit(model_list4, n1)
sem.coefs(model_list4, n1)

# with CN ratios 
model_list5 <- list(
  glm(no_art_5 ~ host + cn_rat + ln.vol_cub, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ host + cn_rat + ln.vol_cub,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + mal + ln.vol_cub, data=n1)
)
sem.fit(model_list5, n1)
sem.coefs(model_list5, n1)

model_list6 <- list(
  glm(no_art_5 ~ host + cn_rat + no_tim_5 + ln.vol_cub, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ host + cn_rat + mal + ln.vol_cub,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + ln.vol_cub, data=n1)
)
sem.fit(model_list6, n1)
sem.coefs(model_list6, n1)

##
model_list7 <- list(
  glm(no_art_5 ~ host + ln.vol_cub + cn_rat, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + host + cn_rat,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + ln.vol_cub + mal, data=n1),
  lm(mal ~ host, data=n1)
)
sem.fit(model_list7, n1)
sem.coefs(model_list7, n1)

model_list8 <- list(
  glm(no_art_5 ~ host + ln.vol_cub + mal, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + host + mal,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + ln.vol_cub + mal, data=n1),
  lm(mal ~ host, data=n1)
)
sem.fit(model_list8, n1)
sem.coefs(model_list8, n1)


model_list9 <- list(
  glm(no_art_5 ~ host + ln.vol_cub + mal, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + host + mal,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + ln.vol_cub, data=n1),
  lm(mal ~ host + cn_rat, data=n1)
)
sem.fit(model_list9, n1)
sem.coefs(model_list9, n1)

model_list10 <- list(
  glm(no_art_5 ~ host + ln.vol_cub + cn_rat, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + host + cn_rat,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + ln.vol_cub, data=n1),
  lm(mal ~ host + no_art_5, data=n1)
)
sem.fit(model_list10, n1)
sem.coefs(model_list10, n1)

model_list11 <- list(
  glm(no_art_5 ~ host + ln.vol_cub + mal, 
      family=quasipoisson, data=n1),
  glm(no_tim_5 ~ ln.vol_cub + host + mal,
      family=quasipoisson, data=n1),
  lm(cn_rat ~ host + ln.vol_cub + no_art_5, data=n1),
  lm(mal ~ host, data=n1)
)
sem.fit(model_list11, n1)
sem.coefs(model_list10, n1)

## CN, mal, abundance path analysis for paper! looks great

# feedback model 1

feed_mod1 <- list(
  
  lm(no_art_5 ~ cn2 + host + ln.vol_cub + con_art5, data=n1), 
  lm(cn2 ~ c.mal + host, data=n1)
  
)

sem.fit(feed_mod1, n1)
coefs12 <- sem.coefs(feed_mod1, n1, standardize="scale")
coefs12
sem.plot(coef.table=coefs12)

# blended feedback model
feed_mod2 <- list(
  
  glm(no_art_5 ~ host + c.mal + ln.vol_cub, data=n1, family=quasipoisson), 
  lm(cn2 ~ host, data=n1),
  lm(c.mal ~ cn2, data=n1),
  lm(mas ~ no_art_5 + ln.vol_cub, data=n1)
  
)

sem.fit(feed_mod2, n1)
coefs2 <- sem.coefs(feed_mod2, n1, standardize="scale")
coefs2
sem.plot(coef.table=coefs12)
##

model_list12 <- list(
  
  lm(mal ~ no_art + host, data=n1),
  lm(cn2 ~ no_art + host + ln.vol_cub, data=n1),
  lm(ln.vol_cub ~ host, data=n1),
  glm(no_art ~ ln.vol_cub + host + c.con_art, family=quasipoisson, data=n1),
  lm(c.con_art ~ host + ln.vol_cub, data=n1),
  lm(slopes ~ no_art + mal, data=n1)
)

sem.fit(model_list12, n1)
coefs12 <- sem.coefs(model_list12, n1, standardize="scale")
coefs12
sem.plot(coef.table=coefs12)

model_list13 <- list(
  
  glm(no_art_5 ~ cn2 + host + ln.vol_cub + c.con_art5, data=n1),
  lm(cn2 ~ mal + host + ln.vol_cub, data=n1),
  lm(ln.vol_cub ~ host, data=n1),
  lm(c.con_art5 ~ host + ln.vol_cub, data=n1)
  
)

sem.fit(model_list13, n1)
coefs13 <- sem.coefs(model_list13, n1, standardize="scale")
sem.plot(coef.table=coefs13)

model_list14 <- list(
  
  glm(no_art ~ cn2 + host + ln.vol_cub, data=n1, family=quasipoisson),
  lm(cn2 ~ mal + host + ln.vol_cub, data=n1),
  lm(ln.vol_cub ~ host, data=n1),
  lm(c.con_art5 ~ host + ln.vol_cub, data=n1), 
  lm(slopes ~ mal + no_art, data=n1)
  
  
)

sem.fit(model_list14, n1)
coefs14 <- sem.coefs(model_list14, n1, standardize="scale")
coefs14
sem.plot(coef.table=coefs14)

model_list15 <- list(
  
  glm(no_art ~ no_tim + host + ln.vol_cub, data=n1, family=quasipoisson),
  glm(no_tim ~ host + ln.vol_cub + con_tim + mal + I(mal^2), data=n1, family=quasipoisson)
 
)

sem.fit(model_list15, n1)
coefs15 <- sem.coefs(model_list15, n1)
coefs15
sem.plot(coef.table=coefs15)

an1 <- lm(slopes ~ no_art_5*host  + mal*host, data=n1[n1$slopes < 0.5 & n1$no_art > 5, ])
summary(an1)

n1.1 <- n1[n1$slopes < 0.5 & n1$no_art > 5, ]
plot(n1.1$cn2, n1.1$slopes)


#### TEST TWO MODELS!

# feedback model
feed_model0 <- list(
  
  glm(no_art_5 ~ cn2 + host + ln.vol_cub, family=quasipoisson, data=n1), 
  lm(cn2 ~ c.mal + host, data=n1), 
  lm(mas ~ cn2 + host + ln.vol_cub + no_art_5, data=n1)
  
)

sem.fit(feed_model0, n1)
coefs_f0 <- sem.coefs(feed_model0, n1)
coefs_f0

# ddep model
ddep_model0 <- list(
  
  glm(no_art_5 ~ ln.vol_cub + host, family=quasipoisson, data=n1), 
  lm(cn2 ~ no_art_5 + host + ln.vol_cub, data=n1), 
  lm(mas ~ host + no_art_5 + ln.vol_cub, data=n1),
  lm(c.mal ~ no_art_5 + host, data=n1)
  
)

sem.fit(ddep_model0, n1)
coefs_d0 <- sem.coefs(ddep_model0, n1)
coefs_d0

# exploratory analysis

exp_model0 <- list(
  
  glm(no_art_5 ~ ln.vol_cub + host, family=quasipoisson, data=n1), 
  lm(cn2 ~ host, data=n1), 
  lm(mas ~ host*c.mal + no_art_5, data=n1),
  lm(c.mal ~ host + cn2, data=n1)
  
)

sem.fit(exp_model0, n1)
coefs_e0 <- sem.coefs(exp_model0, n1)
coefs_e0

exp_model1 <- list(
  
  glm(no_art_5 ~ ln.vol_cub + host, family=quasipoisson, data=n1), 
  lm(cn2 ~ host + c.mal, data=n1), 
  lm(mas ~ host*c.mal + no_art_5, data=n1),
  lm(c.mal ~ host, data=n1)
  
)

sem.fit(exp_model1, n1)
coefs_e1 <- sem.coefs(exp_model1, n1)
coefs_e1

exp_model2 <- list(
  
  glm(no_art_5 ~ ln.vol_cub + host + c.mal, family=quasipoisson, data=n1), 
  lm(cn2 ~ host + c.mal, data=n1), 
  lm(mas ~ host*c.mal + no_art_5, data=n1),
  lm(c.mal ~ host, data=n1)
  
)

sem.fit(exp_model2, n1)
coefs_e2 <- sem.coefs(exp_model2, n1)
coefs_e2
sem.plot(coef.table=coefs_d0)
#

exp_model3 <- list(
  
  glm(no_art_5 ~ ln.vol_cub + host, family=quasipoisson, data=n1), 
  lm(cn2 ~ host + c.mal, data=n1), 
  lm(mas ~ no_art_5 + c.mal, data=n1),
  lm(c.mal ~ host + no_art_5, data=n1)
  
)

sem.fit(exp_model3, n1)
coefs_e3 <- sem.coefs(exp_model3, n1)
coefs_e3
sem.plot(coef.table=coefs_d0)

# test correlated between abundance and abundance - mass slope

am <- mvrnorm(n = 100000, mu=c(0, 0), Sigma=rbind(c(1, .5), c(.5, 1)))
ams <- lapply(1:1000, function(x) mvrnorm(n = 10, mu=c(0, 0), Sigma=rbind(c(1, .5), c(.5, 1))))
ams2 <- do.call(rbind, lapply(ams, function(x) {
  
  an1 <- lm(x[,1] ~ x[,2])
  beta <- coef(an1)[2]
  xbar <- mean(x[,2])
  return(c(beta, xbar))
  
}))

summary(lm(ams2[,1] ~ ams2[, 2]))

# sample from poisson and add

dcoefs <- rep(NA, 100)

for(i in 1:100) {

ma <- rpois(n=100000, lambda = 1)
ma <- ma[ma>0]

counts <- table(ma)
mass <- 1:length(counts)
an1 <- lm(counts ~ mass) # -5470
coef1 <- coef(an1)[2]

probs <- counts / sum(counts)
add <- sample(1:length(counts), size = 10000, prob = probs, replace=TRUE)
counts2 <- table(add)
counts3 <- rep(0, length(counts))
counts3[1:length(counts2)] <- counts2

total <- counts + counts3
an2 <- lm(total ~ mass)
coef2 <- coef(an2)[2]

dcoef <- coef2 - coef1

dcoefs[i] <- dcoef

}

## Using lavaan
n1$sq.no_art_5 <- sqrt(n1$no_art_5)
n1$sq.no_any_5 <- sqrt(n1$no_any_5)

feed_mod1 <- {
'
# regressions
cn2 ~ c.mal * host + c.mal + host
mas ~ sq.no_art_5 
sq.no_art_5 ~ cn2 + ln.vol_cub + con_art5
'
}

fit1 <- sem(feed_mod1, data=n1[n1$host %in% c("A", "C"),])
summary(fit1, standardized=TRUE)

# 
feed_mod2 <- {
  '
  # regressions
  cn2 ~ c.mal + host
  sq.no_any_5 ~ cn2 + ln.vol_cub + con_any5
  '
}

fit2 <- sem(feed_mod2, data=n1[n1$host %in% c("A", "C"),])
summary(fit2, standardized=TRUE)

# feedback with arts and Timema
feed_mod3 <- {
  '
  # regressions
  cn2 ~ c.mal 
  sq.no_any_5 ~ cn2 + ln.vol_cub + con_any5
  '
}

fit3 <- sem(feed_mod3, data=n1[n1$host %in% c("A", "C"),], group="host")
summary(fit3, standardized=TRUE, fit=TRUE, rsquare=TRUE)

# dens dep with arts and Timema
dd_mod1 <- {
  '
  # regressions
  c.mal + cn2 ~ sq.no_any_5
  sq.no_any_5 ~ ln.vol_cub

  # covariancecs
  c.mal ~~ 0*cn2
  '
}

fit4 <- sem(dd_mod1, data=n1[n1$host %in% c("A", "C"),], group="host")
summary(fit4, standardized=TRUE)

# latent variable model for bird predation
feed_mod5 <- {
  '
  # regressions
  bp =~ ln.vol_cub * host
  mas ~ bp + sq.no_art_5 
  sq.no_art_5 ~ bp + ln.vol_cub 
  c.mal ~ bp
  '
}

fit1 <- sem(feed_mod5, data=n1[n1$host %in% c("A", "C"),])
summary(fit1, standardized=TRUE)

## ddep latent variable model
feed_mod6 <- {
  '
  # regressions
  bp =~ no_art_5 + host
  mas ~ bp + no_art_5 
  c.mal ~ bp
  no_art_5 ~ ln.vol_cub
  '
}

fit1 <- sem(feed_mod6, data=n1[n1$host %in% c("A", "C"),])
summary(fit1, standardized=TRUE)

fit2 <- lm(mas ~ no_art_5, data=n1)
summary(fit2)
