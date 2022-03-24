####### 1st December 2014 GMK, APB and TEF #########
####### On the effects of Timema maladaptation on arthropod community structure #####
####### effects of maladaptation on mass-abundance slopes ##########



rm(list=ls())

### calculate mass-abundance relationships for each plant 

dd <- read.csv("/Users/tfarkas/Dropbox/Projects/CA Arthropods/N1/Data/all_art_n1_data.csv")
dd_4 <- subset(dd, dd$is_4==TRUE)
dd_4 <- subset(dd, !is.na(dd$weight) & dd$is_5)
dd_4$plant_id <- as.factor(dd_4$plant_id)
h1 <- hist(dd_4$weight)


# Tim's way -----------------------------------------------------------------------

nd <- data.frame(plant_id=vector(mode="numeric", length=length(levels(dd_4$plant_id))), 
                           mas=vector(mode="numeric", length=length(levels(dd_4$plant_id))))

for (i in 1:length(levels(dd_4$plant_id))) {
  
  hg <- hist(dd_4$weight[dd_4$plant_id==levels(dd_4$plant_id)[i]])
  print(lm(hg$counts ~ hg$mids))
  mod <- lm(hg$counts ~ hg$mids)
  nd$mas[i] <- mod$coefficients[2]
  nd$plant_id[i] <- levels(dd_4$plant_id)[i]
  
}


#### ANDREW's WAY -----------------------------------------------------------------------

# First make the equal spaced bins based on log weight
bins<-20 # you choose!!!
dd_4$mass.bin<-cut(log(dd_4$weight),bins)

# now get the columns that we work with
work<-dd_4[,c("plant_id", "weight", "mass.bin", "plant_sp")]
head(work)
# now loop over each plant ID and do stuff

# plant IDs
pids<-as.numeric(levels(work$plant_id))
ll<-length(pids)

# collection zone
slopes<-rep(NA,ll)
int1 <- rep(NA, ll)
int3 <- rep(NA, ll)
int5 <- rep(NA, ll)
int7 <- rep(NA, ll)
int8 <- rep(NA, ll)
int9 <- rep(NA, ll)
int10 <- rep(NA, ll)
int18 <- rep(NA, ll)
int20 <- rep(NA, ll)
pid <- levels(work$plant_id)
sd <- data.frame(pid, int1, int3, int5, int7, int8, int9, int10, int18, int20, slopes)
for(i in 1:ll){
  # get the plant id subset and check it
  plant<-subset(work, plant_id==pids[i])
  head(plant)
  
  # build the data to be analysed
  # First the data frame with frequencies in each mass bin
  mass.ab<-data.frame(with(plant, table(mass.bin)),cat=levels(plant$mass.bin))
  # force an ascending numberic code on them for plotting
  mass.ab$cat.code0 <- (1:bins) 
  mass.ab$cat.code1 <- mass.ab$cat.code0 - 1 
  mass.ab$cat.code3 <- mass.ab$cat.code0 - 3
  mass.ab$cat.code5 <- mass.ab$cat.code0 - 5
  mass.ab$cat.code7 <- mass.ab$cat.code0 - 7
  mass.ab$cat.code8 <- mass.ab$cat.code0 - 8
  mass.ab$cat.code9 <- mass.ab$cat.code0 - 9
  mass.ab$cat.code10 <- mass.ab$cat.code0 - 10
  mass.ab$cat.code18 <- mass.ab$cat.code0 - 18
  mass.ab$cat.code20 <- mass.ab$cat.code0 - 20
  
  # replace 0 with NAs
  mass.ab$Freq[mass.ab$Freq==0]<-NA
  plot.these<-mass.ab
  
  # # make a plot on the fly (not showing in loop)
  # par(mar=c(8,5,4,4))
  # plot(log(Freq) ~ cat.code, data = plot.these, type="p",
  #      axes=FALSE, xlab="", ylab = "log(Abudance)", ylim=c(0,2))
  # axis(2, at = seq(0,2,0.5))
  # axis(1, at=1:bins, labels=as.character(levels(plant$mass.bin)), las=3)
  
  # run the models
  mod0 <- lm(log(Freq) ~ cat.code0, data = plot.these)
  mod1 <- lm(log(Freq) ~ cat.code1, data = plot.these)
  mod3 <- lm(log(Freq) ~ cat.code3, data = plot.these)
  mod5 <- lm(log(Freq) ~ cat.code5, data = plot.these)
  mod7 <- lm(log(Freq) ~ cat.code7, data = plot.these)
  mod8 <- lm(log(Freq) ~ cat.code8, data = plot.these)
  mod9 <- lm(log(Freq) ~ cat.code9, data = plot.these)
  mod10 <- lm(log(Freq) ~ cat.code10, data = plot.these)
  mod18 <- lm(log(Freq) ~ cat.code18, data = plot.these)
  mod20 <- lm(log(Freq) ~ cat.code20, data = plot.these)
  
  # collect the coefficients
  sd$int1[i] <- coef(mod1)[1]
  sd$int3[i] <- coef(mod3)[1]
  sd$int5[i] <- coef(mod5)[1]
  sd$int7[i] <- coef(mod7)[1]
  sd$int8[i] <- coef(mod8)[1]
  sd$int9[i] <- coef(mod9)[1]
  sd$int10[i] <- coef(mod10)[1]
  sd$int18[i] <- coef(mod18)[1]
  sd$int20[i] <- coef(mod20)[1]
  sd$slopes[i]<-coef(mod0)[2]
  
}

sd$pid <- as.character(sd$pid)
sd$pid <- as.factor(ifelse(nchar(sd$pid) == 1, paste("00", sd$pid, sep = ""), 
                                ifelse(nchar(sd$pid) == 2, paste("0", sd$pid, sep=""), 
                                       sd$pid)))

sd[is.na(sd$slopes), 2:11] <- NA

## add intercepts to n1

n1$mas <- NA
n1$mai1 <- NA
n1$mai3 <- NA
n1$mai5 <- NA
n1$mai7 <- NA
n1$mai8 <- NA
n1$mai9 <- NA
n1$mai10 <- NA
n1$mai18 <- NA
n1$mai20 <- NA


n1[match(sd$pid[sd$pid != "021"], n1$plant_id), 
   c("mai1", "mai3", "mai5", "mai7", "mai8", "mai9", "mai10", "mai18", "mai20", "mas")] <- sd[sd$pid != "021", 2:11]

hist(na.omit(slopes))


slopes

# add slopes to pd

pd <- read.csv("~/Documents/University of Sheffield/PROJECT/Data/pd.csv")

pd$slopes <- rep(NA, 146)


for(i in 1:length(sd$slopes)) {
  for(j in 1:length(pd$plant_id)) {
    if (sd$pid[i]==pd$plant_id[j]) {
      pd$slopes[j] <- sd$slopes[i]; break
    }
  }
}

write.csv(pd, "~/Documents/University of Sheffield/PROJECT/Data/pd.csv")
##### lm mal vs slopes
contrasts(pd$host) <- c(0, 1)

###Adenostoma slopes
anA <- lm(slopes ~ mal + rich_art_4 + no_art_4, data=pd, subset=no_art_4>1&no_tim_4>1&host=="A")
summary(anA)

###Ceanothus slopes
anC <- lm(slopes ~ mal + rich_art_4 + no_art_4, data=pd, subset=no_art_4>1&no_tim_4>1&host=="C")
summary(anC)



plotA <- plot(slopes ~ mal, data=pd, subset=host=="A"&no_art_4>3&no_tim_4>2)
abline(plotA)


### add MAS to plant-level data frame. pd= plant data

pd <- read.csv("~/Documents/University of Sheffield/PROJECT/Data/n1_plant_data_16Sept14.csv")
pd$mas <- rep(NA, 147)
#pd$mas <- vector(length=147, mode="numeric")
pd <- pd[-147,]

for(i in 1:length(nd$plant_id)) {
 for(j in 1:length(pd$plant_id)) {
   
   if (nd$plant_id[i]==pd$plant_id[j]) {
       pd$mas[j] <- nd$mas[i]; break
    }
  }
}

## maladaptation analysis

pd$host <- as.character(pd$host)
pd$host <- as.factor(pd$host)
contrasts(pd$host) <- c(0,1)
an1 <- lm(mas ~ mal*host + ln.vol_cub + no_tim, data=pd)
summary(an1)


# flowering and growth effects 

##new variable. no of flowers added up instead of average for poisson

pd$no_flwrs_total <- pd$no_flwrs_b1 + pd$no_flwrs_b2 + pd$no_flwrs_b3 


hist(pd$avg_growth)
hist(pd$avg_no_flwrs)

an1 <- lm(avg_growth ~ avg_no_flwrs, data=pd)
contrasts(pd$host) <- c(0,1)
summary(an1)

an2 <- glm(no_flwrs_total ~ no_herb_4 + ln.vol_cub + avg_growth, data=pd, family=quasipoisson)
summary(an2)


pd$avg_biomass <- pd$biomass/pd$no_art
pd$avg_biomass

aggregate(art$weight, by=list(art$morph), data=subset(art, plant_sp=="Adenostoma"), FUN=mean, na.rm=TRUE)
names(art)
