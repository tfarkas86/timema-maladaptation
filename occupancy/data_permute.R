rm(list=ls())

dd <- read.csv("C:/Users/tfarkas/Dropbox/Projects/CA Arthropods/AmNatMS/Data/obsTimData2011.csv") # load data
dd$con_mal <- ifelse(dd$host=="A", dd$conC, dd$conA)
dd$con_ad <- ifelse(dd$host=="C", dd$conC, dd$conA)
dd$con_rat <- dd$con_mal/dd$con_ad
dd$pcon_mal <- con_mal/(con)



# patch occupancy w/ con vs. heterospecific connectivity

an1 <- glm(occupied ~ con_mal + con_ad + host + ln.size, family=quasibinomial, data=dd)
summary(an1)

an1 <- glm(occupied ~ conA*host + conC*host + host + ln.size, family=binomial, data=dd)
summary(an1)

# patch occupancy with proportion or ratio heterospecific connectivity

an2 <- glm(occupied ~ con + pcon_mal + host + ln.size, family=quasibinomial, data=dd)
summary(an2) 

an2 <- glm(occupied ~ con_rat + con + host + ln.size, family=quasibinomial, data=dd)
summary(an2) 



