#This script reproduce analysis done in Bartomeus et al. 2017 (submitted)
#On the inconsistency of pollinator species traits for predicting either response to agricultural intensification or functional contribution. 
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
library(mvabund)
library(lattice)

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

#gis range
summary(gis_wt)

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

#gis range
summary(gis_bb)

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

#gis range
summary(gis_cb)

#Response traits----
#Alternative test of fourth corner analysis based in Brown et al and Wang et al papers.
#This is the approach finally used in ms.

#Adapt to wt
# now fit the fourth corner model, only as a function of a few of traits and env variables:
head(gis_wt)
#select only variables for which we have strong predictions
gis_wt <- gis_wt[,c(-3,-4,-7,-8)]

#Option one:
ftSmall=traitglm(R = gis_wt, L = spec_wt, Q = traits_wt[,c(-2,-8,-9,-10,-11)])
ftSmall$fourth.corner
ftSmallcoef <- ftSmall$coefficients
anova(ftSmall, nBoot = 10) #make nBoot=1000 for accutare results
summary(ftSmall) #slow!
#get pseudo-R2 using observed vs predicted
x <- predict.traitglm(object = ftSmall, newR=gis_wt)
plot(unlist(spec_wt) ~ as.vector(x))
summary(lm(unlist(spec_wt) ~ as.vector(x)))

#p.uni
a_wt <- anova(ftSmall, p.uni="adjusted") #make nBoot=1000 for accutare results  #SLOW! 2 h 17 min. gives all p-values.
capture.output(a_wt, file = "data_/a_wt.txt")

#Option 2: with LASSO penealty
ft1=traitglm(R = gis_wt, L = spec_wt, Q = traits_wt[,c(-2,-8,-9,-10,-11)], method="glm1path")
ft1$fourth.corner
ft1coef <- ft1$coefficients
#Check difference between two options
plot(ftSmallcoef, ft1coef[-c(2,55:69,70,75,82,87,94,99,106,108)]) #the shrinkage works also in high coeffs, but presumably with high uncertainty too?
#evaluation of the model
#ft1=traitglm(R = gis_wt[,c(1,4)], L = spec_wt, Q = traits_wt[,c(4,5)], method="glm1path")
x <- predict.traitglm(object = ft1)
plot(unlist(spec_wt) ~ as.vector(x))
summary(lm(unlist(spec_wt) ~ as.vector(x)))

results_wt <- data.frame(crop = "wt", variable = names(ft1$coefficients), estimates = as.numeric(ft1$coefficients))
head(results_wt)

ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
coef(ft1)
plot(ft1)

#plot
#Is neg.bin with log-link, so I need to back-transform.
#As there is standardization of variables, I use standardized values
luc <- scale(gis_wt$ag_300) 
tr <- scale(traits_wt$ITfam) 
resp_q1 <- function(tr, ft1, x){ #units are standardized effect sizes.
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_300")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "ITfam")]*quantile(tr)[2]
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_300:ITfam")]*x*quantile(tr)[2])
}
yq1 <- c()
for(x in c(luc)){
  yq1 <- c(yq1,resp_q1(tr, ft1, x))
}
max(yq1)
resp_q3 <- function(tr, ft1, x){ #standardized effect size
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_300")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "ITfam")]*quantile(tr)[4]
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_300:ITfam")]*x*quantile(tr)[4])
}
yq3 <- c()
for(x in c(luc)){
  yq3 <- c(yq3,resp_q3(tr, ft1, x))
}
max(yq3)
plot(seq(min(c(yq1, yq3)),max(c(yq3, yq1)),length.out = length(luc)) ~ c(luc),
     type = "n", las = 1, xlab = "% Agriculture 300m", ylab = "standardized effect size", xaxt = "n")
axis(side = 1, labels = round(seq(min(gis_wt$ag_300), max(gis_wt$ag_300), length.out = 10)),
     at = seq(min(luc), max(luc), length.out = 10))
temp <- data.frame(luc, yq1)[order(luc),]
lines(x = temp$luc, y = temp$yq1, col = "blue")
temp <- data.frame(luc, yq3)[order(luc),]
lines(x = temp$luc, y = temp$yq3, col = "red")

#plot another interaction
luc <- scale(gis_wt$ag_1500) 
tr <- scale(traits_wt$ITfam) 
resp_q1 <- function(tr, ft1, x){ #units are standardized effect sizes.
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "ITfam")]*quantile(tr)[2]
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500:ITfam")]*x*quantile(tr)[2])
}
yq1 <- c()
for(x in c(luc)){
  yq1 <- c(yq1,resp_q1(tr, ft1, x))
}
max(yq1)
resp_q3 <- function(tr, ft1, x){ #standardized effect size
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "ITfam")]*quantile(tr)[4]
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500:ITfam")]*x*quantile(tr)[4])
}
yq3 <- c()
for(x in c(luc)){
  yq3 <- c(yq3,resp_q3(tr, ft1, x))
}
max(yq3)
plot(seq(min(c(yq1, yq3)),max(c(yq3, yq1)),length.out = length(luc)) ~ c(luc),
     type = "n", las = 1, xlab = "% Agriculture 1500m", ylab = "standardized effect size", xaxt = "n")
axis(side = 1, labels = round(seq(min(gis_wt$ag_1500), max(gis_wt$ag_1500), length.out = 10)),
     at = seq(min(luc), max(luc), length.out = 10))
temp <- data.frame(luc, yq1)[order(luc),]
lines(x = temp$luc, y = temp$yq1, col = "blue")
temp <- data.frame(luc, yq3)[order(luc),]
lines(x = c(luc), y = yq3, col = "red")


#Because all predictors are standardised, you can interpret the size of coefficients as a measure
#of importance. As interaction terms, the fourth coefficients each have an interpretation as the 
#amount by which a unit (1 sd) change in the trait variable changes the slope of the relationship 
#between abundance and a given environmental variable.

#a        = max( abs(ft1$fourth.corner) )
a        = max(0.5)
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)


#BB
spec_bb <- droplevels(spec_bb)
traits_bb <- droplevels(traits_bb)
gis_bb <- gis_bb[,c(-3,-4,-7,-8)]
gis_bb <- droplevels(gis_bb)

ftSmall=traitglm(R = gis_bb, L = spec_bb, Q = traits_bb[,c(-2,-4,-8,-9,-10,-11)])
anova(ftSmall, nBoot = 100) #make nBoot=1000 for accutare results 
#env:trait (fourth corner)    251      64 130.2    0.164 #for 1000 nBoot
a_bb <- anova(ftSmall, p.uni="adjusted") #make nBoot=1000 for accutare results  #SLOW! 2 h 17 min. gives all p-values.
capture.output(a_bb, file = "data_/a_bb.txt")
plot(ftSmall)

ft1=traitglm(R = gis_bb, L = spec_bb, Q = traits_bb[,c(-2,-4, -8,-9,-10,-11)], method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
ft1$coefficients
plot(ft1)

results_bb <- data.frame(crop = "bb", variable = names(ft1$coefficients), estimates = as.numeric(ft1$coefficients))
head(results_bb)

#predictive power
x <- predict.traitglm(object = ft1)
plot(unlist(spec_bb) ~ as.vector(x))
summary(lm(unlist(spec_bb) ~ as.vector(x)))

#a        = max( abs(ft1$fourth.corner) )
a        = max(0.5)
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

#As there is standardization of variables, I use standardized values
luc <- scale(gis_bb$ag_1500) 
tr <- scale(traits_bb$tongue) 
resp_q1 <- function(tr, ft1, x){ #units are standardized effect sizes.
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "tongue")]*quantile(tr)[2]
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500:tongue")]*x*quantile(tr)[2])
}
yq1 <- c()
for(x in c(luc)){
  yq1 <- c(yq1,resp_q1(tr, ft1, x))
}
max(yq1)
resp_q3 <- function(tr, ft1, x){ #standardized effect size
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "tongue")]*quantile(tr)[4]
      + ft1$coefficients[which(names(ft1$coefficients) == "ag_1500:tongue")]*x*quantile(tr)[4])
}
yq3 <- c()
for(x in c(luc)){
  yq3 <- c(yq3,resp_q3(tr, ft1, x))
}
max(yq3)
plot(seq(min(c(yq1, yq3)),max(c(yq3, yq1)),length.out = length(luc)) ~ c(luc),
     type = "n", las = 1, xlab = "% Agriculture 1500m", ylab = "standardized effect size", xaxt = "n")
axis(side = 1, labels = round(seq(min(gis_wt$ag_1500), max(gis_wt$ag_1500), length.out = 10)),
     at = seq(min(luc), max(luc), length.out = 10))
temp <- data.frame(luc, yq1)[order(luc),]
lines(x = temp$luc, y = temp$yq1, col = "blue")
temp <- data.frame(luc, yq3)[order(luc),]
lines(x = temp$luc, y = temp$yq3, col = "red")


#CB
gis_cb <- gis_cb[,c(-3,-4,-7,-8)]

ftSmall=traitglm(R = gis_cb, L = spec_cb, Q = traits_cb[,c(-2,-8,-9,-10,-11)])
anova(ftSmall, nBoot = 100) #make nBoot=1000 for accutare results
a_cb <- anova(ftSmall, p.uni="adjusted") #make nBoot=1000 for accutare results  #SLOW! 2 h 17 min. gives all p-values.
capture.output(a_cb, file = "data_/a_cb.txt")

ft1=traitglm(R = gis_cb, L = spec_cb, Q = traits_cb[,c(-2,-8,-9,-10,-11)], method="glm1path")
ft1$coefficients
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
plot(ft1)
#predictive power
x <- predict.traitglm(object = ft1)
plot(unlist(spec_cb) ~ as.vector(x))
summary(lm(unlist(spec_cb) ~ as.vector(x)))

results_cb <- data.frame(crop = "cb", variable = names(ft1$coefficients), estimates = as.numeric(ft1$coefficients))
head(results_cb)


#a        = max( abs(ft1$fourth.corner) )
a        = max(0.5)
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

#plot interaction
luc <- scale(gis_cb$open_300) 
tr <- scale(traits_cb$PDrar20) 
resp_q1 <- function(tr, ft1, x){ #units are standardized effect sizes.
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "open_300")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "PDrar20")]*quantile(tr)[2]
      + ft1$coefficients[which(names(ft1$coefficients) == "open_300:PDrar20")]*x*quantile(tr)[2])
}
yq1 <- c()
for(x in c(luc)){
  yq1 <- c(yq1,resp_q1(tr, ft1, x))
}
max(yq1)
resp_q3 <- function(tr, ft1, x){ #standardized effect size
  exp(ft1$coefficients[which(names(ft1$coefficients) == "")] 
      + ft1$coefficients[which(names(ft1$coefficients) == "open_300")]*x
      + ft1$coefficients[which(names(ft1$coefficients) == "PDrar20")]*quantile(tr)[4]
      + ft1$coefficients[which(names(ft1$coefficients) == "open_300:PDrar20")]*x*quantile(tr)[4])
}
yq3 <- c()
for(x in c(luc)){
  yq3 <- c(yq3,resp_q3(tr, ft1, x))
}
max(yq3)
plot(seq(min(c(yq1, yq3)),max(c(yq3, yq1)),length.out = length(luc)) ~ c(luc),
     type = "n", las = 1, xlab = "% semi-natural 300m", ylab = "standardized effect size", xaxt = "n")
axis(side = 1, labels = round(seq(min(gis_wt$ag_1500), max(gis_wt$ag_1500), length.out = 10)),
     at = seq(min(luc), max(luc), length.out = 10))
temp <- data.frame(luc, yq1)[order(luc),]
lines(x = temp$luc, y = temp$yq1, col = "blue")
temp <- data.frame(luc, yq3)[order(luc),]
lines(x = c(luc), y = yq3, col = "red")

#unify results
results <- rbind(results_wt, results_bb, results_cb)
write.csv(results, file = "data_/coefs.csv")


#fourth corner analysis----
#this analysis is not used in the final paper as using mvabund is a better approach
#However, you can see here results are consistent.

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
#plot(four6b_bb, stat = "D2", alpha = 0.1)

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


#Are species abundance responding to LUC?----
#One of reviewers ask, so we check. (not in final paper)

gis_wt2 <- gis_wt
gis_wt2$abund <- rowSums(spec_wt)
head(gis_wt2)
m <- lm(abund ~ ag_300, data = gis_wt2)
#plot(m)
summary(m)
m <- lm(abund ~ ag_1500, data = gis_wt2)
summary(m)
m <- lm(abund ~ open_300, data = gis_wt2)
summary(m) #0.09
#test
plot(predict(m) ~ gis_wt2$abund)
summary(lm(predict(m) ~ gis_wt2$abund)) #interesting...

m <- lm(abund ~ open_1500, data = gis_wt2)
summary(m) #0.06
plot(gis_wt2$abund ~ gis_wt2$open_1500)

#bb
gis_bb$abund <- rowSums(spec_bb)
head(gis_bb)
m <- lm(abund ~ ag_300, data = gis_bb)
#plot(m)
summary(m) #0.04
m <- lm(abund ~ ag_1500, data = gis_bb)
summary(m)
m <- lm(abund ~ open_300, data = gis_bb)
summary(m) 
m <- lm(abund ~ open_1500, data = gis_bb)
summary(m) 
plot(gis_bb$abund ~ gis_bb$ag_300)
abline(m$coefficients + 0.30)

#cb
gis_cb$abund <- rowSums(spec_cb)
head(gis_cb)
m <- lm(abund ~ ag_300, data = gis_cb)
#plot(m)
summary(m)
m <- lm(abund ~ ag_1500, data = gis_cb)
summary(m) #0.01
m <- lm(abund ~ open_300, data = gis_cb)
summary(m) 
m <- lm(abund ~ open_1500, data = gis_cb)
summary(m) #0.06
