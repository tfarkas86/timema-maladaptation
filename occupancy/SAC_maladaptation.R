rm(list=ls())

dd <- read.csv("PATH/obsTimData2011.csv") # load data

# does alternate-host connectivity (pcon_mal) influence the likelihood that a 
# Timema is maladapted?

an1 <- lm(mal ~ host + pcon_mal + total + ln.size + con, data=dd)
summary(an1) # seems like yes, but the pattern goes opposite to expectations...
plot(an1) # ick, but the residuals show strong patterns. something is missing

# should probably use binomial anyway

b_mal <- cbind(dd$no_mal, dd$no_cam) # binomal response variable for maladaptation

an2 <- glmer(b_mal ~ pcon_mal + host + total + ln.size + con + (1|id), family=binomial, data=dd) 
summary(an2) # seems like no here, but host is still significant

# but what about spatial autocorrelation? here is the GLMM PQL method, except
# the model is already mixed, needing id as a random factor

dd_long <- cnt2bin(dd, "no_mal", "no_cam") # use my awesome function to make data long form

an3 <- glmmPQL(bin ~ pcon_mal + host, random=~1|id, 
               correlation=corExp(form=~x+y),
               family=binomial, data=dd_long)
summary(an3)

# ugh. since many individual Timema have the same location (are on the same bush)
# this correlation function does not work. we could fake it by adding jitter to one 
# of the coordinates

an4 <- glmmPQL(bin ~ pcon_mal + host, random=~1|id, 
               correlation=corExp(form=~jitter(x)+y),
               family=binomial, data=dd_long)
summary(an4) # and it runs okay, showing no big difference from the last model

# but is the jitter method okay? should something else be done? im not convinced
# that this is an appropriate way to get around this problem ...

