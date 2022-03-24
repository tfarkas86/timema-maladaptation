#### create distance matrix for all bushes using rdist in field package
rm(list=ls())

n1 <- read.csv("C:/Users/tfarkas/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_plant_data_16Sept14.csv")
n1 <- n1[-147, ]
loc.mat <- cbind(n1$x, n1$y)
library(fields)
dist.mat <- rdist(loc.mat)
disp <- 2 # dispersal distance is 2m 
alpha <- 1/disp

##### baby problem #####

loc.mat.b <- loc.mat[1:3,]
dist.mat.b <- rdist(loc.mat.b)
diag(dist.mat.b) <- 0
con.mat.b <- (alpha^2/2*pi) * exp(-1 * alpha *dist.mat.b)
diag(con.mat.b) <- 0
host.b <- c("A", "A", "B")
con.2N.b <- rowSums(con.mat.b[,host.b=="A"])

#######

is.nan.data.frame <- function(x) # function to remove NaNs
  do.call(cbind, lapply(x, is.nan))

con <- data.frame(n1$plant_id) # initialize a data frame to output (with host IDs)

con.mat.2N <- (alpha^2/2*pi) * exp(-1 * alpha *dist.mat) # connectivity matrix, no weights
diag(con.mat.2N) <- 0 
con.mat.2N <- signif(con.mat.2N, 5)

con$con_2N <- rowSums(con.mat.2N) # connectivity, no population size weights
con$con_2N_A <- rowSums(con.mat.2N[,n1$host=="A"]) # unweighted connectivity to Adenostoma
con$con_2N_C <- rowSums(con.mat.2N[,n1$host=="C"]) # unweighted connectivity to Ceanothus

n1$volume <- log(4/3 * pi * (n1$height/2) * (n1$width/2) * (n1$length/2)) # elliptical volume
con.mat.vol <- t(t(con.mat.2N)*n1$volume) # connectivity matrix weighted by plant volume
con$con_vol <- rowSums(con.mat.vol) # volume-weighted connectivity
con$con_vol_A <- rowSums(con.mat.vol[,n1$host=="A"]) # VWC Adenostoma
con$con_vol_C <- rowSums(con.mat.vol[,n1$host=="C"]) # VWC Ceanothus

# timema only
con.mat.tim <- t(t(con.mat.2N)*n1$no_tim) # Timema popualtion-weighted connectivity matrix
con$con_tim <- rowSums(con.mat.tim)
con$con_tim_A <- rowSums(con.mat.tim[,n1$host=="A"])
con$con_tim_C <- rowSums(con.mat.tim[,n1$host=="C"])

# phenotype dependent : qimmA and qimmC et al.

con.mat.gp <- t(t(con.mat.2N)*n1$no_tim * (n1$no_grn_tim/n1$no_tim_gs)) # gp = proprotion green timema
con.mat.gp[is.nan(con.mat.gp)] <- NA # remove NaNs
con$con_gp_A <- rowSums(con.mat.gp[,n1$host=="A"], na.rm=TRUE) # intermediate to qimmA
con$con_gp_C <- rowSums(con.mat.gp[,n1$host=="C"], na.rm=TRUE) # ................qimmC


con$q2N_A <- con$con_gp_A/con$con_tim_A # qimmA: expected green frequency from A
con$q2N_C <- con$con_gp_C/con$con_tim_C # qimmC: ..............................C

con.mat.grn <- t(t(con.mat.2N)*n1$no_grn_tim) # weighted matrix by green timema only
con.mat.str <- t(t(con.mat.2N)*n1$no_str_tim) # .................  striped .........


con$con_grn_A <- rowSums(con.mat.grn[,n1$host=="A"], na.rm=TRUE)  # connectivity weighted by green timema on A
con$con_grn_C <- rowSums(con.mat.grn[,n1$host=="C"], na.rm=TRUE)  # .........................................C
con$con_str_A <- rowSums(con.mat.str[,n1$host=="A"], na.rm=TRUE)  # .........................striped.........A
con$con_str_C <- rowSums(con.mat.str[,n1$host=="C"], na.rm=TRUE)  # . .......................................C



# other Timema coinnectivities

con.mat.tim4 <- t(t(con.mat.2N)*n1$no_tim_4) 
con.mat.tim5 <- t(t(con.mat.2N)*n1$no_tim_5) 
con.mat.tim6 <- t(t(con.mat.2N)*n1$no_tim_6) 
con.mat.tim7 <- t(t(con.mat.2N)*n1$no_tim_7) 
con.mat.tim8 <- t(t(con.mat.2N)*n1$no_tim_8) 
con.mat.tim9 <- t(t(con.mat.2N)*n1$no_tim_9) 
con.mat.tim10 <- t(t(con.mat.2N)*n1$no_tim_10) 
con$con_tim4 <- rowSums(con.mat.tim4)
con$con_tim5 <- rowSums(con.mat.tim5)
con$con_tim6 <- rowSums(con.mat.tim6)
con$con_tim7 <- rowSums(con.mat.tim7)
con$con_tim8 <- rowSums(con.mat.tim8)
con$con_tim9 <- rowSums(con.mat.tim9)
con$con_tim10 <- rowSums(con.mat.tim10)

# all arthropods, including timema

con.mat.any <- t(t(con.mat.2N)*n1$no_any) 
con.mat.any4 <- t(t(con.mat.2N)*n1$no_any_4) 
con.mat.any5 <- t(t(con.mat.2N)*n1$no_any_5)
con.mat.any6 <- t(t(con.mat.2N)*n1$no_any_6) 
con.mat.any7 <- t(t(con.mat.2N)*n1$no_any_7) 
con.mat.any8 <- t(t(con.mat.2N)*n1$no_any_8) 
con.mat.any9 <- t(t(con.mat.2N)*n1$no_any_9) 
con.mat.any10 <- t(t(con.mat.2N)*n1$no_any_10) 
con$con_any <- rowSums(con.mat.any)
con$con_any4 <- rowSums(con.mat.any4)
con$con_any5 <- rowSums(con.mat.any5)
con$con_any6 <- rowSums(con.mat.any6)
con$con_any7 <- rowSums(con.mat.any7)
con$con_any8 <- rowSums(con.mat.any8)
con$con_any9 <- rowSums(con.mat.any9)
con$con_any10 <- rowSums(con.mat.any10)


# non-Timema arthropods
con.mat.art <- t(t(con.mat.2N)*n1$no_art) 
con.mat.art4 <- t(t(con.mat.2N)*n1$no_art_4) 
con.mat.art5 <- t(t(con.mat.2N)*n1$no_art_5) 
con.mat.art6 <- t(t(con.mat.2N)*n1$no_art_6) 
con.mat.art7 <- t(t(con.mat.2N)*n1$no_art_7) 
con.mat.art8 <- t(t(con.mat.2N)*n1$no_art_8) 
con.mat.art9 <- t(t(con.mat.2N)*n1$no_art_9) 
con.mat.art10 <- t(t(con.mat.2N)*n1$no_art_10) 
con$con_art <- rowSums(con.mat.art)
con$con_art4 <- rowSums(con.mat.art4)
con$con_art5 <- rowSums(con.mat.art5)
con$con_art6 <- rowSums(con.mat.art6)
con$con_art7 <- rowSums(con.mat.art7)
con$con_art8 <- rowSums(con.mat.art8)
con$con_art9 <- rowSums(con.mat.art9)
con$con_art10 <- rowSums(con.mat.art10)

# herbivores, including Timema
con.mat.herb <- t(t(con.mat.2N)*n1$no_herb) 
con.mat.herb4 <- t(t(con.mat.2N)*n1$no_herb_4) 
con.mat.herb5 <- t(t(con.mat.2N)*n1$no_herb_5) 
con.mat.herb6 <- t(t(con.mat.2N)*n1$no_herb_6) 
con.mat.herb7 <- t(t(con.mat.2N)*n1$no_herb_7) 
con.mat.herb8 <- t(t(con.mat.2N)*n1$no_herb_8) 
con.mat.herb9 <- t(t(con.mat.2N)*n1$no_herb_9) 
con.mat.herb10 <- t(t(con.mat.2N)*n1$no_herb_10) 
con$con_herb <- rowSums(con.mat.herb)
con$con_herb4 <- rowSums(con.mat.herb4)
con$con_herb5 <- rowSums(con.mat.herb5)
con$con_herb6 <- rowSums(con.mat.herb6)
con$con_herb7 <- rowSums(con.mat.herb7)
con$con_herb8 <- rowSums(con.mat.herb8)
con$con_herb9 <- rowSums(con.mat.herb9)
con$con_herb10 <- rowSums(con.mat.herb10)

# non-Timema herbivores
con.mat.art.herb <- t(t(con.mat.2N)*n1$no_art_herb) 
con.mat.art.herb4 <- t(t(con.mat.2N)*n1$no_art_herb_4) 
con.mat.art.herb5 <- t(t(con.mat.2N)*n1$no_art_herb_5) 
con.mat.art.herb6 <- t(t(con.mat.2N)*n1$no_art_herb_6) 
con.mat.art.herb7 <- t(t(con.mat.2N)*n1$no_art_herb_7) 
con.mat.art.herb8 <- t(t(con.mat.2N)*n1$no_art_herb_8) 
con.mat.art.herb9 <- t(t(con.mat.2N)*n1$no_art_herb_9) 
con.mat.art.herb10 <- t(t(con.mat.2N)*n1$no_art_herb_10) 
con$con_art_herb <- rowSums(con.mat.art.herb)
con$con_art_herb4 <- rowSums(con.mat.art.herb4)
con$con_art_herb5 <- rowSums(con.mat.art.herb5)
con$con_art_herb6 <- rowSums(con.mat.art.herb6)
con$con_art_herb7 <- rowSums(con.mat.art.herb7)
con$con_art_herb8 <- rowSums(con.mat.art.herb8)
con$con_art_herb9 <- rowSums(con.mat.art.herb9)
con$con_art_herb10 <- rowSums(con.mat.art.herb10)

# carnivores
con.mat.carn <- t(t(con.mat.2N)*n1$no_carn) 
con.mat.carn4 <- t(t(con.mat.2N)*n1$no_carn_4) 
con.mat.carn5 <- t(t(con.mat.2N)*n1$no_carn_5) 
con.mat.carn6 <- t(t(con.mat.2N)*n1$no_carn_6) 
con.mat.carn7 <- t(t(con.mat.2N)*n1$no_carn_7) 
con.mat.carn8 <- t(t(con.mat.2N)*n1$no_carn_8) 
con.mat.carn9 <- t(t(con.mat.2N)*n1$no_carn_9) 
con.mat.carn10 <- t(t(con.mat.2N)*n1$no_carn_10) 
con$con_carn <- rowSums(con.mat.carn)
con$con_carn4 <- rowSums(con.mat.carn4)
con$con_carn5 <- rowSums(con.mat.carn5)
con$con_carn6 <- rowSums(con.mat.carn6)
con$con_carn7 <- rowSums(con.mat.carn7)
con$con_carn8 <- rowSums(con.mat.carn8)
con$con_carn9 <- rowSums(con.mat.carn9)
con$con_carn10 <- rowSums(con.mat.carn10)

# carnivores
con.mat.lep <- t(t(con.mat.2N)*n1$no_lep) 
con.mat.lep4 <- t(t(con.mat.2N)*n1$no_lep_4) 
con.mat.lep5 <- t(t(con.mat.2N)*n1$no_lep_5) 
con.mat.lep6 <- t(t(con.mat.2N)*n1$no_lep_6) 
con.mat.lep7 <- t(t(con.mat.2N)*n1$no_lep_7) 
con.mat.lep8 <- t(t(con.mat.2N)*n1$no_lep_8) 
con.mat.lep9 <- t(t(con.mat.2N)*n1$no_lep_9) 
con.mat.lep10 <- t(t(con.mat.2N)*n1$no_lep_10) 
con$con_lep <- rowSums(con.mat.lep)
con$con_lep4 <- rowSums(con.mat.lep4)
con$con_lep5 <- rowSums(con.mat.lep5)
con$con_lep6 <- rowSums(con.mat.lep6)
con$con_lep7 <- rowSums(con.mat.lep7)
con$con_lep8 <- rowSums(con.mat.lep8)
con$con_lep9 <- rowSums(con.mat.lep9)
con$con_lep10 <- rowSums(con.mat.lep10)

write.csv(con, file="C:/Users/tfarkas/Dropbox/Projects/CA Arthropods/N1/Data/csv_for_r/n1_connectivities.csv")
