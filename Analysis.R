#This script reproduce analysis done in Bartomeus et al. 2016 (submitted)
#Pollinator species traits do not predict either response to agricultural 
# intensification or functional contribution 
#Ignasi Bartomeus, Daniel P. Cariveau, Tina Harrison, Rachael Winfree.

#load data----

net <- read.csv("data_/specimens.csv")
traits <- read.csv("data_/traits.csv")
gis <- read.csv("data_/gis.csv")


#load libraries----

library(reshape)
library(ade4)
library(FD)
library(MuMIn)

#perpare data for each crop----

##wt 
#spec
net_wt <- subset(net, Crop == "wt")
net_wt <- droplevels(net_wt)
#str(net_wt)
spec_wt <- cast(net_wt, Farm ~ Gen_sp, fun = length, value = "Round")
head(spec_wt)
rownames(spec_wt) <- spec_wt$Farm
spec_wt <- spec_wt[,-1]
#gis
gis_wt <- subset(gis, Crop == "wt", select = -c(Crop))
rownames(gis_wt) <- gis_wt$Farm
gis_wt <- subset(gis_wt, select=-c(Farm))
#traits
traits_wt <- subset(traits, Gen_sp %in% colnames(spec_wt))
str(traits_wt)
rownames(traits_wt) <- traits_wt$Gen_sp
traits_wt <- traits_wt[,c(-1)]

##bb 
#spec
net_bb <- subset(net, Crop == "bb")
net_bb <- droplevels(net_bb)
str(net_bb)
spec_bb <- cast(net_bb, Farm ~ Gen_sp, fun = length, value = "Round")
head(spec_bb)
rownames(spec_bb) <- spec_bb$Farm
spec_bb <- spec_bb[,-1]
#gis 
gis_bb <- subset(gis, Crop=="bb", select=-c(Crop))
rownames(gis_bb) <- gis_bb$Farm
gis_bb <- subset(gis_bb, select=-c(Farm))
str(gis_bb)
#traits
traits_bb <- subset(traits, Gen_sp %in% colnames(spec_bb))
str(traits_bb)
rownames(traits_bb) <- traits_bb$Gen_sp
traits_bb <- traits_bb[,c(-1)]

##cb
#spec
net_cb <- subset(net, Crop == "cb")
net_bb <- droplevels(net_cb)
str(net_cb)
spec_cb <- cast(net_cb, Farm ~ Gen_sp, fun = length, value = "Round")
head(spec_cb)
rownames(spec_cb) <- spec_cb$Farm
spec_cb <- spec_cb[,-1]
str(spec_cb)
#gis 
gis_cb <- subset(gis, Crop=="cb", select=-c(Crop))
rownames(gis_cb) <- gis_cb$Farm
gis_cb <- subset(gis_cb, select=-c(Farm))
str(gis_cb)
#traits
traits_cb <- subset(traits, Gen_sp %in% colnames(spec_cb))
str(traits_cb)
rownames(traits_cb) <- traits_cb$Gen_sp
traits_cb <- traits_cb[,c(-1)]

#fourth corner analysis----

##wt
four6b_wt <- fourthcorner(tabR = gis_wt ,tabL = spec_wt ,tabQ = traits_wt[,c(-2,-8,-9,-10,-11)],
                          nrepet = 999, modeltype = 6, p.adjust.method.G = 'none', p.adjust.method.D = 'none')
print(four6b_wt, stat = "D2")
plot(four6b_wt, stat = "D2", alpha = 0.05)

#visualization in Fig 1
#calculate CWM ~ Ag300
source("functcomp_fixed.R")
temp <- as.matrix(spec_wt)
colnames(temp) <- colnames(spec_wt)
rownames(temp) <- rownames(spec_wt)
fc <- functcomp_fixed(traits_wt[,c(-2,-8,-9,-10,-11)], a = temp)
temp <- cbind(gis_wt, fc)
plot(temp$ag_300, temp$ITfam, las = 1, xlab = "% agriculture at 300m buffer", ylab = "body size")
abline(summary(lm(temp$ITfam ~ temp$ag_300)))

##cb
four6b_cb <- fourthcorner(tabR = gis_cb ,tabL = spec_cb ,tabQ = traits_cb[,c(-2,-8,-9,-10,-11)],
                          nrepet = 999, modeltype = 6, p.adjust.method.G = 'none', p.adjust.method.D = 'none')
print(four6b_cb, stat = "D2")
plot(four6b_cb, stat = "D2", alpha = 0.05)
#plot(four6b_cb, stat = "D2", alpha = 0.1)

##bb
four6b_bb <- fourthcorner(tabR = gis_bb,tabL = spec_bb ,tabQ = traits_bb[,c(-2,-8,-9,-10,-11)],
                          nrepet = 999, modeltype = 6, p.adjust.method.G = 'none', p.adjust.method.D = 'none')
print(four6b_bb, stat = "D2")
plot(four6b_bb, stat = "D2", alpha = 0.05)

#efficiency analysis----

#1) Visitation:----
head(traits_wt)
m <- lm(visitation_wt ~ Nest_place + Sociality + Parasitic + ITfam + PDrar20 + tongue, 
        data = traits_wt, na.action = na.fail)
summary(m)
out <- dredge(m)
subset(out, delta <2)
importance(out)
model.avg(out, subset= delta < 2, revised.var = TRUE) 

m <- lm(visitation_wt ~ tongue, data = traits_wt, na.action = na.fail)
summary(m)
#plot(m)
plot(traits_wt$visitation_wt ~ traits_wt$tongue)
abline(m) #totally driven by bombus impatiens!

#cb
m <- lm(visitation_cb ~ Nest_place + Sociality + Parasitic + ITfam + PDrar20 + tongue, 
        data = traits_cb, na.action = na.fail)
summary(m)
out <- dredge(m)
subset(out, delta <2)
importance(out)
model.avg(out, subset= delta < 2, revised.var = TRUE) 

#best is nest only
m <- lm(visitation_cb ~ Sociality + tongue, data = traits_cb, na.action = na.fail)
summary(m)
m <- lm(visitation_cb ~ Nest_place, data = traits_cb, na.action = na.fail)
summary(m)
anova(m)
#plot(m)
plot(traits_cb$visitation_cb ~ traits_cb$Nest_place, las = 1, ylab = "visitation rate", xlab = "nest site")
plot(traits_cb$visitation_cb ~ traits_cb$Sociality)
plot(traits_cb$visitation_cb ~ traits_cb$tongue)

#bb
m <- lm(visitation_bb ~ Nest_place + Sociality + ITfam + PDrar20 + tongue, 
        data = traits_bb, na.action = na.fail)
summary(m)
out <- dredge(m)
subset(out, delta <2)
importance(out)
#model.avg(out, subset= delta < 2, revised.var = TRUE) 

#best is PDrar20 only
m <- lm(visitation_bb ~ PDrar20, data = traits_bb, na.action = na.fail)
summary(m)
#plot(m)
plot(traits_bb$visitation_bb ~ traits_bb$PDrar20, las = 1, xlab = "diet specialism", 
     ylab = "visitation rate")
abline(m) #driven by Andrena specialists. Still interesting.

#2)Efficiency.----

#load group level data
eff_wt3  <-  read.csv("data_/efficiency_wt.csv", header=TRUE)

m <- lm(pol_mean ~ Nest_place + Sociality + ITfam + PDrar20 + tongue, data = eff_wt3, na.action = na.fail)
summary(m)
out <- dredge(m)
subset(out, delta <2)
importance(out)

#best is ITfam & tongue sig only
m <- lm(pol_mean ~ ITfam + tongue, data = eff_wt3, na.action = na.fail)
summary(m)
m <- lm(pol_mean ~ ITfam, data = eff_wt3, na.action = na.fail)
summary(m)
#without bombus
#eff_wt4 <- eff_wt3[-1,]
#m <- lm(pol_mean ~ ITfam, data = eff_wt4, na.action = na.fail)
#summary(m)

#plot(m)
plot(eff_wt3$pol_mean ~ eff_wt3$ITfam, xlab = "body size", ylab = "pollination efficiency", las = 1)
abline(m) 
#plot(eff_wt4$pol_mean ~ eff_wt4$ITfam)
#abline(m) #not driven by bombus

#cb
#load group level data
eff_cb3  <-  read.csv("data_/efficiency_cb.csv", header=TRUE)

m <- lm(pol_mean ~ Nest_place + Sociality + ITfam + PDrar20 + tongue, data = eff_cb3, na.action = na.fail)
summary(m)
out <- dredge(m)
subset(out, delta <2)
importance(out)

#best is null! tongue next, let's plot it
m <- lm(pol_mean ~ ., data = eff_cb3, na.action = na.fail)
summary(m)
m <- lm(pol_mean ~ tongue, data = eff_cb3, na.action = na.fail)
summary(m)
#plot(m)
plot(eff_cb3$pol_mean ~ eff_cb3$tongue)
abline(m) #as expected, no pattern!

#bb
#load group level data
eff_bb3  <-  read.csv("data_/efficiency_bb.csv", header=TRUE)

m <- lm(pol_mean ~ ITfam + PDrar20 + tongue, data = eff_bb3, na.action = na.fail) #nest place and soc removed due to low variability
summary(m)
out <- dredge(m)
subset(out, delta <2)
importance(out)

m <- lm(pol_mean ~ ITfam + PDrar20 + tongue, data = eff_bb3, na.action = na.fail)
summary(m)
#but watch out:
cor(eff_bb3$ITfam, eff_bb3$tongue)
m <- lm(pol_mean ~ tongue, data = eff_bb3, na.action = na.fail)
summary(m)
plot(eff_bb3$pol_mean ~ eff_bb3$tongue, xlab = "tongue length", ylab = "pollination efficiency", las = 1)
abline(m)
m <- lm(pol_mean ~ ITfam , data = eff_bb3, na.action = na.fail)
plot(eff_bb3$pol_mean ~ eff_bb3$ITfam)
abline(m)
#plot(m)



