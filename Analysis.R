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

#Sup Mat: Alternative test of fourth corner analysis based in Brown et al and Wang et al papers.----

#Adapt to wt
# now fit the fourth corner model, only as a function of a couple of traits and env variables:
ftSmall=traitglm(R = gis_wt, L = spec_wt, Q = traits_wt[,c(-2,-8,-9,-10,-11)])
anova(ftSmall, nBoot = 10) #make nBoot=1000 for accutare results
anova(ftSmall, p.uni="adjusted") #make nBoot=1000 for accutare results  #SLOW! 2 h 17 min. gives all p-values.
# Model: traitglm(L = spec_wt, R = gis_wt, Q = traits_wt[, c(-2, -8, -9, 
#                                                            Model:     -10, -11)])
# 
# Multivariate test:
#   Res.Df Df.diff   Dev Pr(>Dev)   
# (Intercept)                             764                          
# V1                                      764       0  0.00     1.00   
# sppAgapostemon_texanus                  763       1  6.04     0.60   
# sppAgapostemon_virescens                762       1  4.71     0.61   
# sppAugochlora_pura                      761       1 13.63     0.50   
# sppAugochlorella_aurata                 760       1  0.00     0.93   
# sppAugochloropsis_metallica             759       1  9.00     0.57   
# sppBombus_bimaculatus                   758       1  9.27     0.61   
# sppBombus_griseocollis                  757       1  9.56     0.65   
# sppBombus_impatiens                     756       1 59.81     0.57   
# sppCalliopsis_andreniformis             755       1  6.17     0.63   
# sppCeratina_calcarata_dupla_miqmaki     754       1 10.35     0.49   
# sppCeratina_strenua                     753       1  0.00     0.97   
# sppHalictus_confusus                    752       1  0.00     0.91   
# sppHalictus_ligatus                     751       1  0.38     0.63   
# sppHalictus_parallelus                  750       1  9.59     0.62   
# sppHalictus_rubicundus                  749       1  3.30     0.54   
# sppLasioglossum_admirandum              748       1 10.16     0.52   
# sppLasioglossum_albipenne               747       1 10.56     0.59   
# sppLasioglossum_bruneri                 746       1  5.28     0.53   
# sppLasioglossum_coriaceum               745       1 11.31     0.60   
# sppLasioglossum_cressonii               744       1  9.56     0.47   
# sppLasioglossum_ephialtum               743       1  6.20     0.54   
# sppLasioglossum_illinoense              742       1  3.01     0.58   
# sppLasioglossum_imitatum                741       1 60.42     0.45   
# sppLasioglossum_leucocomum              740       1  2.87     0.57   
# sppLasioglossum_mitchelli               739       1  3.94     0.54   
# sppLasioglossum_nymphaearum             738       1  4.08     0.63   
# sppLasioglossum_oblongum                737       1 12.10     0.63   
# sppLasioglossum_obscurum                736       1 10.17     0.61   
# sppLasioglossum_paradmirandum           735       1  0.01     0.92   
# sppLasioglossum_pectorale               734       1  2.68     0.64   
# sppLasioglossum_pilosum                 733       1  6.59     0.51   
# sppLasioglossum_rozeni                  732       1 13.23     0.61   
# sppLasioglossum_smilacinae              731       1 14.22     0.63   
# sppLasioglossum_tegulare                730       1  2.44     0.52   
# sppLasioglossum_trigeminum              729       1 10.94     0.68   
# sppLasioglossum_truncatum               728       1 10.37     0.60   
# sppLasioglossum_versatum                727       1 18.93     0.58   
# sppLasioglossum_weemsi                  726       1  0.12     0.68   
# sppLasioglossum_zephyrum                725       1  6.08     0.52   
# sppMegachile_mendica                    724       1 14.18     0.47   
# sppMelissodes_bimaculata                723       1  3.06     0.67   
# sppPeponapis_pruinosa                   722       1  5.05     0.48   
# sppTriepeolus_remigatus                 721       1 26.52     0.43   
# sppXylocopa_virginica                   720       1  0.16     0.82   
# ag_300                                  719       1  2.33     0.42   
# open_300                                718       1 14.23     0.02 * 
#   forest_edge_300                         717       1  0.45     0.58   
# shannonH_300                            716       1  1.18     0.48   
# ag_1500                                 715       1  0.85     0.56   
# open_1500                               714       1  2.64     0.24   
# forest_edge_1500                        713       1  0.07     0.84   
# shannonH_1500                           712       1  0.01     0.96   
# ag_300.squ                              711       1  0.44     0.63   
# open_300.squ                            710       1  3.88     0.18   
# forest_edge_300.squ                     709       1  0.92     0.49   
# shannonH_300.squ                        708       1  0.46     0.64   
# ag_1500.squ                             707       1  2.81     0.25   
# open_1500.squ                           706       1  5.63     0.07 . 
# forest_edge_1500.squ                    705       1  0.43     0.56   
# shannonH_1500.squ                       704       1 17.36     0.01 **
#   ag_300.Nest_placehole                   703       1  0.01     0.76   
# ag_300.Nest_placesoil                   702       1  0.33     0.65   
# ag_300.Nest_placestem                   701       1  4.93     0.04 * 
#   ag_300.Nest_placewood                   700       1  0.86     0.27   
# ag_300.Socialityfac_social              699       1  0.19     0.71   
# ag_300.SocialitySolitary                698       1  3.65     0.17   
# ag_300.ParasiticYes                     697       1  4.11     0.14   
# ag_300.ITfam                            696       1  3.37     0.10 . 
# ag_300.PDrar20                          695       1  1.06     0.54   
# ag_300.tongue                           694       1  0.04     0.90   
# open_300.Nest_placehole                 693       1  3.64     0.11   
# open_300.Nest_placesoil                 692       1  2.22     0.30   
# open_300.Nest_placestem                 691       1  0.28     0.40   
# open_300.Nest_placewood                 690       1  0.00     0.97   
# open_300.Socialityfac_social            689       1  0.00     0.93   
# open_300.SocialitySolitary              688       1  0.28     0.72   
# open_300.ParasiticYes                   687       1  0.98     0.46   
# open_300.ITfam                          686       1  3.56     0.05 * 
#   open_300.PDrar20                        685       1  1.38     0.44   
# open_300.tongue                         684       1  0.08     0.81   
# forest_edge_300.Nest_placehole          683       1  2.54     0.07 . 
# forest_edge_300.Nest_placesoil          682       1  3.68     0.18   
# forest_edge_300.Nest_placestem          681       1  0.08     0.70   
# forest_edge_300.Nest_placewood          680       1  0.07     0.68   
# forest_edge_300.Socialityfac_social     679       1  0.05     0.78   
# forest_edge_300.SocialitySolitary       678       1  1.78     0.33   
# forest_edge_300.ParasiticYes            677       1  1.97     0.32   
# forest_edge_300.ITfam                   676       1  0.17     0.59   
# forest_edge_300.PDrar20                 675       1  0.58     0.62   
# forest_edge_300.tongue                  674       1  2.02     0.20   
# shannonH_300.Nest_placehole             673       1  0.02     0.15   
# shannonH_300.Nest_placesoil             672       1  1.21     0.36   
# shannonH_300.Nest_placestem             671       1  1.50     0.19   
# shannonH_300.Nest_placewood             670       1  0.11     0.67   
# shannonH_300.Socialityfac_social        669       1  1.04     0.25   
# shannonH_300.SocialitySolitary          668       1  4.74     0.09 . 
# shannonH_300.ParasiticYes               667       1  1.60     0.27   
# shannonH_300.ITfam                      666       1  4.36     0.08 . 
# shannonH_300.PDrar20                    665       1  0.11     0.78   
# shannonH_300.tongue                     664       1  2.87     0.12   
# ag_1500.Nest_placehole                  663       1  0.00     0.35   
# ag_1500.Nest_placesoil                  662       1  0.00     0.98   
# ag_1500.Nest_placestem                  661       1  0.27     0.56   
# ag_1500.Nest_placewood                  660       1  0.39     0.51   
# ag_1500.Socialityfac_social             659       1  0.97     0.34   
# ag_1500.SocialitySolitary               658       1  0.28     0.65   
# ag_1500.ParasiticYes                    657       1  0.00     1.00   
# ag_1500.ITfam                           656       1  3.76     0.05 * 
#   ag_1500.PDrar20                         655       1  0.00     0.99   
# ag_1500.tongue                          654       1 53.94     0.02 * 
#   open_1500.Nest_placehole                653       1  0.00     0.38   
# open_1500.Nest_placesoil                652       1  0.00     0.99   
# open_1500.Nest_placestem                651       1 55.88     0.02 * 
#   open_1500.Nest_placewood                650       1  0.57     0.38   
# open_1500.Socialityfac_social           649       1  0.14     0.71   
# open_1500.SocialitySolitary             648       1  8.52     0.06 . 
# open_1500.ParasiticYes                  647       1  2.95     0.10 . 
# open_1500.ITfam                         646       1  0.17     0.70   
# open_1500.PDrar20                       645       1  4.93     0.13   
# open_1500.tongue                        644       1  2.48     0.20   
# forest_edge_1500.Nest_placehole         643       1  0.00     0.44   
# forest_edge_1500.Nest_placesoil         642       1  0.02     0.94   
# forest_edge_1500.Nest_placestem         641       1  0.00     1.00   
# forest_edge_1500.Nest_placewood         640       1 53.63     0.03 * 
#   forest_edge_1500.Socialityfac_social    639       1  0.00     0.99   
# forest_edge_1500.SocialitySolitary      638       1 -0.06     1.00   
# forest_edge_1500.ParasiticYes           637       1 69.98     0.02 * 
#   forest_edge_1500.ITfam                  636       1  0.48     0.41   
# forest_edge_1500.PDrar20                635       1  0.00     1.00   
# forest_edge_1500.tongue                 634       1  1.44     0.29   
# shannonH_1500.Nest_placehole            633       1  0.01     0.03 * 
#   shannonH_1500.Nest_placesoil            632       1 55.08     0.02 * 
#   shannonH_1500.Nest_placestem            631       1  2.17     0.15   
# shannonH_1500.Nest_placewood            630       1  0.27     0.50   
# shannonH_1500.Socialityfac_social       629       1  1.37     0.29   
# shannonH_1500.SocialitySolitary         628       1  0.06     0.78   
# shannonH_1500.ParasiticYes              627       1  0.00     0.79   
# shannonH_1500.ITfam                     626       1  1.55     0.28   
# shannonH_1500.PDrar20                   625       1  0.94     0.55   
# shannonH_1500.tongue                    624       1  0.61     0.42   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Arguments: P-value calculated using 99 resampling iterations via PIT-trap block resampling (to account for correlation in testing).



ft1=traitglm(R = gis_wt, L = spec_wt, Q = traits_wt[,c(-2,-8,-9,-10,-11)], method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
plot(ft1)
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
gis_bb <- droplevels(gis_bb)

ftSmall=traitglm(R = gis_bb, L = spec_bb, Q = traits_bb[,c(-2,-4,-8,-9,-10,-11)])
anova(ftSmall, nBoot = 100) #make nBoot=1000 for accutare results 
#env:trait (fourth corner)    251      64 130.2    0.164 #for 1000 nBoot


plot(ftSmall)

ft1=traitglm(R = gis_bb, L = spec_bb, Q = traits_bb[,c(-2,-4, -8,-9,-10,-11)], method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
plot(ft1)

#a        = max( abs(ft1$fourth.corner) )
a        = max(0.5)
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

#CB
ftSmall=traitglm(R = gis_cb, L = spec_cb, Q = traits_cb[,c(-2,-8,-9,-10,-11)])
anova(ftSmall, nBoot = 100) #make nBoot=1000 for accutare results

ft1=traitglm(R = gis_cb, L = spec_cb, Q = traits_cb[,c(-2,-8,-9,-10,-11)], method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
anova(ft1, nBoot = 100) #make nBoot=1000 for accutare results
plot(ft1)

#without shannon/edges
ftSmall=traitglm(R = gis_cb[,c(-3,-4,-7,-8)], L = spec_cb, Q = traits_cb[,c(-2,-8,-9,-10,-11)])
anova(ftSmall, nBoot = 1000) #make nBoot=1000 for accutare results
#sig!!
ft1=traitglm(R = gis_cb[,c(-3,-4,-7,-8)], L = spec_cb, Q = traits_cb[,c(-2,-8,-9,-10,-11)], method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
plot(ft1)


#a        = max( abs(ft1$fourth.corner) )
a        = max(0.5)
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

#Are species abundance responding to LUC?----

gis_wt$abund <- rowSums(spec_wt)
head(gis_wt)
m <- lm(abund ~ ag_300, data = gis_wt)
#plot(m)
summary(m)
m <- lm(abund ~ ag_1500, data = gis_wt)
summary(m)
m <- lm(abund ~ open_300, data = gis_wt)
summary(m) #0.09
m <- lm(abund ~ open_1500, data = gis_wt)
summary(m) #0.06

#Looks like it does
