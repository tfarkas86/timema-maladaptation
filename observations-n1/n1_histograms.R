rm(list=ls())

# import data
n1 <- read.csv("~/Dropbox/Projects/CA Arthropods/Data/2013/Network1/n1_data.csv")

# make histogram composite

pdf(file="./n1_histograms.pdf", width=10, height=10)
par(mfrow=c(2,2))

# timema histogram Adenostoma
tim_hist <- hist(n1$no_tim[n1$plant_sp=="Adenostoma"], breaks=seq(from=min(n1$no_tim)-.5, to=max(n1$no_tim)+.5, by=1), 
                 xlim=c(-.5, 20.5), ylim=c(0, 50), xlab="no. Timema", ylab="no. bushes", 
                 main="Timema histogram, Adenostoma")

# artropod histogram Adenostoma

art_hist <- hist(n1$no_art_4[n1$plant_sp=="Adenostoma"], breaks=seq(from=min(n1$no_art_4)-.5, to=max(n1$no_art_4)+.5, by=1), 
                 xlim=c(-.5, 60.5), ylim=c(0, 20), xlab="no. arts > 4mm", ylab="no. bushes", 
                 main="arthropod histogram, Adenostoma")

# timema histogram Ceanothus
tim_hist <- hist(n1$no_tim[n1$plant_sp=="Ceanothus"], breaks=seq(from=min(n1$no_tim)-.5, to=max(n1$no_tim)+.5, by=1), 
                 xlim=c(-.5, 20.5), ylim=c(0, 50), xlab="no. Timema", ylab="no. bushes", 
                 main="Timema histogram, Ceanothus")

# artropod histogram Ceanothus

art_hist <- hist(n1$no_art_4[n1$plant_sp=="Ceanothus"], breaks=seq(from=min(n1$no_art_4)-.5, to=max(n1$no_art_4)+.5, by=1), 
                 xlim=c(-.5, 60.5), ylim=c(0, 20), xlab="no. arts > 4mm", ylab="no. bushes", 
                 main="arthropod histogram, Ceanothus")
dev.off()

