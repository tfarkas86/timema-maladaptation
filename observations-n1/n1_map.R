rm(list=ls())

n1 <- read.csv("~/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data_16Sept14.csv")
n1 <- n1[-147,]
n1$host <- factor(n1$host)
n1$c.con_art5 <- n1$con_art5-mean(n1$con_art5)
n1$c.mal <- n1$mal-mean(n1$mal, na.rm=TRUE)
n1$occ <- ifelse(n1$no_tim_gs==0, FALSE, TRUE)
n1$radius <- (n1$vol_cub/(100000*pi))^(1/3)

n1_occ <- n1[n1$occ,]

adOnly <- n1_occ[n1_occ$host=="A",]
ceOnly <- n1_occ[n1_occ$host=="C",]
empty <- n1[!n1$occ,]

cFunA <- colorRamp(c("skyblue1", "blue"))
cFunC <- colorRamp(c("darkorange2", "peachpuff"))

colsA <- rgb(cFunA(adOnly$mal), maxColorValue=255)
colsC <- rgb(cFunC(ceOnly$mal), maxColorValue=255)

dev.off()
par(fg="lightgrey")
symbols(x=adOnly$x, y=adOnly$y, circles=(adOnly$radius), bg=colsA, inches=FALSE, xlim=c(0,75), ylim=c(-1, 50),
        xlab = "metres E-W", ylab="metres N-S", fg="darkgreen", lwd=1.5, cex.lab=1.5, cex.axis=1.5)
symbols(x=ceOnly$x, y=ceOnly$y, circles=(ceOnly$radius), bg=colsC, inches=F, add=T, fg="darkred", lwd=1.5)
symbols(x=empty$x, y=empty$y, circles=(empty$radius), bg="white", fg=ifelse(empty$host=="A", "darkgreen", "darkred"),
        inches=F, add=T, lwd=1.5)


# # test for relationship
# 
# data <- data.frame(pointID=obsDataAll$PlantID, x=obsDataAll$X, y=obsDataAll$Y,
#                    area = obsDataAll$area, abun=obsDataAll$total, occ=as.logical(x=obsDataAll$occupied))
# data <- transform(data, occ=ifelse(data$abun==0, FALSE, TRUE))
# 
# # contruct pixelating function 
# pixT <- 50 # total number of pixels
# xMax <- max(data$x)
# xMin <- min(data$x)
# yMax <- max(data$y)
# yMin <- min(data$y)
# xRange<- diff(range(data$x))
# yRange <- diff(range(data$y))
# if(xRange > yRange) ratio <- xRange/yRange else
#   ratio <- yRange/xRange
# pixX <- round(sqrt(pixT/ratio)) # approximate # pixels for x axis
# 
# # find proper multiple for number of pixels
# i = 0
# repeat {
#   rem1 <- pixT%%(pixX+i)
#   rem2 <- pixT%%(pixX-i)
#   if(rem1 == 0) {
#     pixX <- pixX+i 
#     print(pixX)
#     break}
#   if(rem2 == 0) {
#     pixX <- pixX-i 
#     print(pixX)
#     break}
#   i <- i+1
# }
# 
# pixY <- pixT/pixX
# 
# # define boundaries for pixels
# 
# xLims <- seq(from=xMin, to=xMax, length.out=pixX+1)
# yLims <- seq(from=yMin, to=yMax, length.out=pixY+1)
# 
# # draw a grid
# i=0
# for (i in 1:(pixX+1)) {
#   lines(x=rep(xLims[i], 2), y=c(min(yLims), max(yLims)))}
# for (i in 1:(pixY+1)) {
#   lines(y=rep(yLims[i], 2), x=c(min(xLims), max(xLims)))}
# 
# # assign sites to pixels
# 
# pixDataXmin <- rep(xLims[-(pixX+1)], pixY)
# pixDataXmax <- rep(xLims[-1], pixY)
# 
# pixDataYmin <- NULL
# pixDataYmax <- NULL
# for (i in 1:(pixY)) {
#   pixDataYmin <- c(pixDataYmin, rep(yLims[i], pixX))
#   pixDataYmax <- c(pixDataYmax, rep(yLims[i+1], pixX))
# }
# 
# pixData <- data.frame(pixel=1:pixT, xMin=pixDataXmin, xMax=pixDataXmax,
#                       yMin=pixDataYmin, yMax=pixDataYmax)
# pointData <- data
# 
# # loop to check inclusion of data points
# for(i in 1:length(pointData$pointID)) { for(j in 1:pixT) { 
#   if(pointData[i,"x"] >= pixData[j, "xMin"] &&
#        pointData[i,"x"] <= pixData[j, "xMax"] &&
#        pointData[i,"y"] >= pixData[j, "yMin"] &&
#        pointData[i,"y"] <= pixData[j, "yMax"]) {
#     
#     pointData[i, "pixel"] <- pixData[j, "pixel"]   
#   }}}
# 
# # data <- data.frame(pointID=obsDataAll$PlantID, x=obsDataAll$X, y=obsDataAll$Y,
# #                    area = obsDataAll$area, abun=obsDataAll$total, occ=as.logical(x=obsDataAll$occupied))
# # data <- transform(data, occ=ifelse(data$abun==0, FALSE, TRUE))
# # 
# # noCells <- 25
# # max <- max(c(data$x, data$y))
# # min <- min(c(data$x, data$y))
# # n <- max/sqrt(noCells)
# # sqLims <- seq(from=min, to=max, n)
# # 
# # # draw a grid
# # 
# # #dev.off(); points(data$x, data$y, xlim=c(-1,max), ylim=c(-1, max))
# # for (i in 1:length(sqLims)) {
# #   lines(x=rep(sqLims[i], 2), y=c(min(sqLims), max(sqLims)))
# #   lines(y=rep(sqLims[i], 2), x=c(min(sqLims), max(sqLims))) }
# # 
# # # which points in which cells?
# # pixDataXmin <- rep(sqLims[-length(sqLims)], sqrt(noCells))
# # pixDataXmax <- rep(sqLims[-1], sqrt(noCells))
# # 
# # pixDataYmin <- vector()
# # pixDataYmax <- vector()
# # for (i in 1:(length(sqLims)-1)) {
# #   pixDataYmin <- c(pixDataYmin, rep(sqLims[i], sqrt(noCells)))
# #   pixDataYmax <- c(pixDataYmax, rep(sqLims[i+1], sqrt(noCells)))
# # }
# # 
# # pixData <- data.frame(pixel=1:noCells, xMin=pixDataXmin, xMax=pixDataXmax,
# #                       yMin=pixDataYmin, yMax=pixDataYmax)
# # dataOut <- data.frame(data)
# # 
# # for(i in 1:length(data$pointID)) { for(j in 1:noCells) { # loop to check inclusion of data points
# #   if(data[i,"x"] >= pixData[j, "xMin"] &&
# #        data[i,"x"] <= pixData[j, "xMax"] &&
# #        data[i,"y"] >= pixData[j, "yMin"] &&
# #        data[i,"y"] <= pixData[j, "yMax"]) {
# #     
# #     dataOut[i, "pixel"] <- pixData[j, "pixel"]   
# #   }}}
# # 
# # dataOut <- dataOut[order(dataOut$pixel),] # sort data by pixel number
# 
# pixData$totAbun <- rep(0, noCells)
# pixData$occ <- rep(0, noCells)
# pixData$unOcc <- rep(0, noCells)
# pixData$occArea <- rep(0, noCells)
# pixData$unOccArea <- rep(0, noCells)
# 
# for (i in 1:noCells) { for (j in 1:length(dataOut$pointID)) {
#   if(pixData[i, "pixel"] == dataOut[j, "pixel"]) {
#     pixData[i, "totAbun"] <- pixData[i, "totAbun"] + dataOut[j, "abun"]    
#     if(dataOut[j,"occ"]) {
#       pixData[i, "occ"] <- pixData[i, "occ"] + 1
#       pixData[i, "occArea"] <- pixData[i, "occArea"] + dataOut[j,"area"]} else {
#         pixData[i, "unOcc"] <- pixData[i, "unOcc"] + 1
#         pixData[i, "unOccArea"] <- pixData[i, "unOccArea"] + dataOut[j, "area"]}
#   }}} 
# 
# pixData$noSites <- pixData$occ + pixData$unOcc
# pixData$propOcc <- pixData$occ/(pixData$noSites)
# pixData$avgAbun <- pixData$totAbun/pixData$occ
# pixData$dens <- pixData$totAbun/pixData$occArea
# 
# pixData <- pixData[pixData$noSites > 2,]
# 
# # pixData <- pixData[pixData$noSites<30,]
# 
# 
# an1 <- lm(propOcc ~ dens, data=pixData)
# plot(propOcc ~ dens, data=pixData)
# abline(an1)
# summary(an1)
# an2 <- glm(cbind(occ, unOcc) ~ dens, data=pixData, family=binomial)
# summary(an2)
# 
# data <- pixelAssign(data=data, n=50, plot=T)
# 
# ######### GREEN PHENOTYPE BUBBLE PLOT ############
# 
# obsData1 <- read.table("obsTimDataOcc1.txt", header=T)
# obsData1$area <- (.5*obsData1$length/100)*(.5*obsData1$width/100)*pi
# obsData1$rad <- sqrt(obsData1$area/pi)
# 
# cFun <- colorRamp(c("lightgreen", "darkgreen"))
# colsAll <- rgb(cFun(obsData1$greenProNRG), maxColorValue=255)
# 
# symbols(x=obsData1$X, y=obsData1$Y, circles=(obsData1$rad), bg=colsAll, cex.axis=1.2,
#         xlab="metres", ylab="metres", inches=F, xlim=c(0,35), ylim=c(-1, 60))
# 
# ######## OCCUPIED BUBBLE PLOT ############
# rm(list=ls())
# setwd("~/Dropbox/Projects/CA Arthropods/Data for R")
# obsDataAll <- read.table("obsTimDataAllBushes.txt", header=T)
# obsData1 <- read.table("obsTimDataOcc1.txt", header=T)
# 
# obsDataAll$area <- (.5*obsDataAll$Length/100)*(.5*obsDataAll$Width/100)*pi
# obsDataAll$rad <- sqrt(obsDataAll$area/pi)
# 
# 
# 
# 
# colour <- ifelse(obsDataAll$occupied==0, "white", "black")
# 
# 
# par(bg="white")
# symbols(x=obsDataAll$X, y=obsDataAll$Y, circles=(obsDataAll$rad), bg=colour, inches=F, xlim=c(0,35), ylim=c(-1, 60),
#         xlab = "metres", ylab="metres", cex.lab=1.7, cex.axis=1.7, lwd=2)
# legend(x=20, y=62, legend=c("occupied", "unoccupied"), pch=c(16, 21), col=c("black", "black"), 
#        bty="n", cex=1.8, y.intersp=.5, pt.lwd=2)
# 
# ####### HOST BUBBLE PLOT #######
# 
# rm(list=ls())
# setwd("~/Dropbox/Projects/CA Arthropods/Data for R")
# obsDataAll <- read.table("obsTimDataAllBushes.txt", header=T)
# 
# obsDataAll$area <- (.5*obsDataAll$Length/100)*(.5*obsDataAll$Width/100)*pi
# obsDataAll$rad <- sqrt(obsDataAll$area/pi)
# colour <- ifelse(obsDataAll$host=="A","blue", "orange")
# 
# par(mar=c(rep(4, 4)))
# symbols(x=obsDataAll$X, y=obsDataAll$Y, circles=(obsDataAll$rad), bg=colour, inches=F, xlim=c(0,35), ylim=c(-1, 60),
#         xlab = "metres", ylab="metres", cex.lab=1.2, cex.axis=1.2)
# legend(x=15, y=62, legend=c("Adenostoma", "Ceanothus"), pch=c(16), col=c("blue", "orange"), 
#        bty="n", cex=.9, y.intersp=.8)


