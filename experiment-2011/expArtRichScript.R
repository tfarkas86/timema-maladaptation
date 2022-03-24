# load and attach data

richData <- read.table("artData2011.txt", header = T)
attach(data)

# create contrast codes for striped vs. green comparison

stVgr <- rep(c(0,0,-.5,.5), 5)
x2 <- rep(c(0,1,-.5,-.5), 5)
x3 <- rep(c(1.5, -.5, -.5, -.5), 5)

an1 <- lm(dens5Sp ~ stVgr + x2 + x3)
summary(an1)

# create contrast codes for striped vs. Ceonothus

stVCeo <- rep(c(1,0,0,-1),5)
x2 <- rep(c(-1,0,2,-1),5)
x3 <- rep(c(-1,3,-1,-1),5)

an1 <- lm(dens5Sp ~ stVCeo + x2 + x3)
summary(an1)

# create contrast codes for well vs. poorly adapted

wVp <- rep(c(1,0,1,-2),5)
x2 <- rep(c(-1,0,1,0),5)
x3 <- rep(c(-1,3,-1,-1),5)

an1 <- lm(dens5Sp ~ wVp + x2 + x3)
summary(an1)

ceoRich <- c(richData$dens5Sp[1], richData$dens5Sp[5], richData$dens5Sp[9], 
            richData$dens5Sp[13], richData$dens5Sp[17])
conRich <- c(richData$dens5Sp[2], richData$dens5Sp[6], richData$dens5Sp[10], 
            richData$dens5Sp[14], richData$dens5Sp[18])
grnRich <- c(richData$dens5Sp[3], richData$dens5Sp[7], richData$dens5Sp[11], 
            richData$dens5Sp[15], richData$dens5Sp[19])
strRich <- c(richData$dens5Sp[4], richData$dens5Sp[8], richData$dens5Sp[12], 
            richData$dens5Sp[16], richData$dens5Sp[20])

artRichOut <- data.frame("control"=conRich, "Ceonothus"=ceoRich, "green"=grnRich, "striped"=strRich)
write.table(artRichOut, "artRichPrism.txt")
read.table("artRichPrism.txt", header=T)