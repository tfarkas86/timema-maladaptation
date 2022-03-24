####### 1st December 2014 GMK, APB and TEF #########
####### On the effects of Timema maladaptation on arthropod community structure #####
####### effects of maladaptation on mass-abundance slopes ##########



rm(list=ls())

### calculate mass-abundance relationships for each plant 

dd <- read.csv("/Users/tfarkas/Dropbox/Projects/CA Arthropods/N1/Data/all_art_n1_data.csv", 
               colClasses=c(plant_id = "factor"))
mad <- droplevels(mad[!is.na(mad$weight) & mad$is_4 & mad$plant_id != "021", ])

mas <- aggregate(dd$weight, by=list(dd$plant_id), FUN=function(x) {
  
  if(length(x) < 3) return(NA)
  lw <- log(x)
  h1 <- hist(x)
  cnt <- h1$counts
  cnt[cnt==0] <- NA
  mids <- h1$mids
  an1 <- lm(log(cnt) ~ log(mids))
  return(an1$coefficients[2])
  
})


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
work<-dd_4[,c(1,2,8,70)]
head(work)



# now loop over each plant ID and do stuff

# plant IDs
pids<-as.numeric(levels(work$plant_id))
ll<-length(pids)

# collection zone
slopes<-rep(NA,ll)
pid <- levels(work$plant_id)
sd <- data.frame(slopes, pid)
for(i in 1:ll){
  # get the plant id subset and check it
  plant<-subset(work, plant_id==pids[i])
  head(plant)
  
  # build the data to be analysed
  # First the data frame with frequencies in each mass bin
  mass.ab<-data.frame(with(plant, table(mass.bin)),cat=levels(plant$mass.bin))
  # force an ascending numberic code on them for plotting
  mass.ab$cat.code<-1:bins
  # replace 0 with NAs
  mass.ab$Freq[mass.ab$Freq==0]<-NA
  plot.these<-mass.ab
  
  # make a plot on the fly (not showing in loop)
  par(mar=c(8,5,4,4))
  plot(log(Freq) ~ cat.code, data = plot.these, type="p",
       axes=FALSE, xlab="", ylab = "log(Abudance)", ylim=c(0,2))
  axis(2, at = seq(0,2,0.5))
  axis(1, at=1:bins, labels=as.character(levels(plant$mass.bin)), las=3)
  
  # run the model
  mod<-lm(log(Freq) ~ cat.code, data = plot.these)
  
  # collect the coefficient
  sd$slopes[i]<-coef(mod)[2]
  
}

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
