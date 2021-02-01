#-----
# KB Gorman
# Pygoscelis Frontiers Manuscript
# Isotopic niche analysis
# 31 Jan, 2021
#-----

rm(list=ls())
graphics.off()

setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

install.packages("devtools")
library(devtools)
devtools::install_github("andrewljackson/SIBER@v2.1.5", build_vingettes = TRUE)
library(SIBER)

ALL.Peng.DataSet<- read.table("ALL_Peng_2019-3.txt",header=T,sep="\t") # This is the one that has all
# the penguin data as each community by species, year, and location.
ALL.Peng.DataSet

my.siber.data<- createSiberObject(ALL.Peng.DataSet)
my.siber.data

#Running the bayesian form of the model
# options for running jags
parms <- list()
parms$n.iter <- 10^6   # number of iterations to run the model for, 1,000,000
parms$n.burnin <- 5*10^4 # discard the first set of values, 50,000
parms$n.thin <- 50     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

install.packages("rjags")
library(rjags)

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.

ellipses.posterior <- siberMVN(my.siber.data, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group", 
                 prn=TRUE)

#Output from 1/31/21
Probability values for Column 1 
Mode 0.0434 Mean 0.051 Median 0.0486 
95 % lower = 0.0273 upper = 0.0789 
75 % lower = 0.0333 upper = 0.0618 
50 % lower = 0.0375 upper = 0.0538 
Probability values for Column 2 
Mode 0.124 Mean 0.151 Median 0.141 
95 % lower = 0.0687 upper = 0.254 
75 % lower = 0.0861 upper = 0.186 
50 % lower = 0.1 upper = 0.157 
Probability values for Column 3 
Mode 0.23 Mean 0.267 Median 0.255 
95 % lower = 0.143 upper = 0.415 
75 % lower = 0.175 upper = 0.324 
50 % lower = 0.196 upper = 0.282 
Probability values for Column 4 
Mode 0.147 Mean 0.167 Median 0.16 
95 % lower = 0.093 upper = 0.256 
75 % lower = 0.111 upper = 0.201 
50 % lower = 0.124 upper = 0.176 
Probability values for Column 5 
Mode 0.0395 Mean 0.0439 Median 0.042 
95 % lower = 0.0242 upper = 0.0672 
75 % lower = 0.0292 upper = 0.053 
50 % lower = 0.0327 upper = 0.0463 
Probability values for Column 6 
Mode 0.316 Mean 0.35 Median 0.335 
95 % lower = 0.193 upper = 0.534 
75 % lower = 0.233 upper = 0.423 
50 % lower = 0.261 upper = 0.37 
Probability values for Column 7 
Mode 0.116 Mean 0.131 Median 0.125 
95 % lower = 0.0717 upper = 0.201 
75 % lower = 0.0873 upper = 0.158 
50 % lower = 0.0976 upper = 0.139 
Probability values for Column 8 
Mode 0.0496 Mean 0.0572 Median 0.0545 
95 % lower = 0.0306 upper = 0.089 
75 % lower = 0.0374 upper = 0.0691 
50 % lower = 0.0421 upper = 0.0604 
Probability values for Column 9 
Mode 0.0579 Mean 0.0651 Median 0.062 
95 % lower = 0.0342 upper = 0.102 
75 % lower = 0.0415 upper = 0.0789 
50 % lower = 0.0473 upper = 0.0687 
Probability values for Column 10 
Mode 0.241 Mean 0.252 Median 0.245 
95 % lower = 0.158 upper = 0.357 
75 % lower = 0.185 upper = 0.297 
50 % lower = 0.201 upper = 0.267 
Probability values for Column 11 
Mode 0.224 Mean 0.245 Median 0.238 
95 % lower = 0.151 upper = 0.353 
75 % lower = 0.177 upper = 0.29 
50 % lower = 0.195 upper = 0.26 
Probability values for Column 12 
Mode 0.175 Mean 0.192 Median 0.186 
95 % lower = 0.119 upper = 0.275 
75 % lower = 0.139 upper = 0.227 
50 % lower = 0.153 upper = 0.203 
Probability values for Column 13 
Mode 0.109 Mean 0.124 Median 0.118 
95 % lower = 0.0662 upper = 0.192 
75 % lower = 0.0814 upper = 0.15 
50 % lower = 0.0911 upper = 0.131 

#-----
#Graphing the SEAb

install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)

DS1<- subset(ALL.FW.DataSet, community==1)
DS2<- subset(ALL.FW.DataSet, community==2)
DS3<- subset(ALL.FW.DataSet, community==3)
DS4<- subset(ALL.FW.DataSet, community==4)
DS5<- subset(ALL.FW.DataSet, community==5)
DS6<- subset(ALL.FW.DataSet, community==6)
DS7<- subset(ALL.FW.DataSet, community==7)
DS8<- subset(ALL.FW.DataSet, community==8)
DS9<- subset(ALL.FW.DataSet, community==9)
DS10<- subset(ALL.FW.DataSet, community==10)
DS11<- subset(ALL.FW.DataSet, community==11)
DS12<- subset(ALL.FW.DataSet, community==12)
DS13<- subset(ALL.FW.DataSet, community==13)

my.siber.data1<- createSiberObject(DS1)
my.siber.data2<- createSiberObject(DS2)
my.siber.data3<- createSiberObject(DS3)
my.siber.data4<- createSiberObject(DS4)
my.siber.data5<- createSiberObject(DS5)
my.siber.data6<- createSiberObject(DS6)
my.siber.data7<- createSiberObject(DS7)
my.siber.data8<- createSiberObject(DS8)
my.siber.data9<- createSiberObject(DS9)
my.siber.data10<- createSiberObject(DS10)
my.siber.data11<- createSiberObject(DS11)
my.siber.data12<- createSiberObject(DS12)
my.siber.data13<- createSiberObject(DS13)


# Start here with graphing 

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

par(mfrow=c(4,3), mar=c(5.1,4.1,4.1,2.1))

plot1SO<- plotSiberObject(my.siber.data1,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'))

my.siber.data
c.id<- 1
g.id<- 1

plot1<- plot(iso2~iso1, data=DS1, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="black", main="2008")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

plot2SO<- plotSiberObject(my.siber.data2,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'))

my.siber.data
c.id<- 2
g.id<- 1

plot2<- plot(iso2~iso1, data=DS2, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="darkorange", main="2008")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkorange",
                     lty = 1,
                     lwd = 2)

plot3SO<- plotSiberObject(my.siber.data3,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'),)

c.id<- 3
g.id<- 1

plot3<- plot(iso2~iso1, data=DS3, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="darkblue", main="2008")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkblue",
                     lty = 1,
                     lwd = 2)

plot4SO<- plotSiberObject(my.siber.data4,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'),)

c.id<- 4
g.id<- 1

plot4<- plot(iso2~iso1, data=DS4, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="black", main="2009")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

plot5SO<- plotSiberObject(my.siber.data5,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'),)

c.id<- 5
g.id<- 1

plot5<- plot(iso2~iso1, data=DS5, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="darkorange", main="2009")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkorange",
                     lty = 1,
                     lwd = 2)

plot6SO<- plotSiberObject(my.siber.data6,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'),)

c.id<- 6
g.id<- 1

plot6<- plot(iso2~iso1, data=DS6, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="darkblue", main="2009")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkblue",
                     lty = 1,
                     lwd = 2)

plot7SO<- plotSiberObject(my.siber.data7,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'),)

c.id<- 7
g.id<- 1

plot7<- plot(iso2~iso1, data=DS7, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="black", main="2010")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

plot8SO<- plotSiberObject(my.siber.data8,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'))

c.id<- 8
g.id<- 1

plot8<- plot(iso2~iso1, data=DS8, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="darkorange", main="2010")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkorange",
                     lty = 1,
                     lwd = 2)

plot9SO<- plotSiberObject(my.siber.data9,
                          ax.pad = 2, 
                          hulls = F, community.hulls.args, 
                          ellipses = T, group.ellipses.args,
                          group.hulls = F, group.hull.args,
                          bty = "L",
                          iso.order = c(1,2),
                          xlab = expression({delta}^13*C~'\u2030'),
                          ylab = expression({delta}^15*N~'\u2030'))

c.id<- 9
g.id<- 1

plot9<- plot(iso2~iso1, data=DS9, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
             pch=16, col="darkblue", main="2010")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkblue",
                     lty = 1,
                     lwd = 2)

plot10SO<- plotSiberObject(my.siber.data10,
                           ax.pad = 2, 
                           hulls = F, community.hulls.args, 
                           ellipses = T, group.ellipses.args,
                           group.hulls = F, group.hull.args,
                           bty = "L",
                           iso.order = c(1,2),
                           xlab = expression({delta}^13*C~'\u2030'),
                           ylab = expression({delta}^15*N~'\u2030'))

c.id<- 10
g.id<- 1

plot10<- plot(iso2~iso1, data=DS10, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
              pch=16, col="black", main="2008")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

plot11SO<- plotSiberObject(my.siber.data11,
                           ax.pad = 2, 
                           hulls = F, community.hulls.args, 
                           ellipses = T, group.ellipses.args,
                           group.hulls = F, group.hull.args,
                           bty = "L",
                           iso.order = c(1,2),
                           xlab = expression({delta}^13*C~'\u2030'),
                           ylab = expression({delta}^15*N~'\u2030'))

c.id<- 11
g.id<- 1

plot11<- plot(iso2~iso1, data=DS11, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
              pch=16, col="black", main="2009")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

plotSO12<- plotSiberObject(my.siber.data12,
                           ax.pad = 2, 
                           hulls = F, community.hulls.args, 
                           ellipses = T, group.ellipses.args,
                           group.hulls = F, group.hull.args,
                           bty = "L",
                           iso.order = c(2,1),
                           xlab = expression({delta}^13*C~'\u2030'),
                           ylab = expression({delta}^15*N~'\u2030'))

c.id<- 12
g.id<- 1

plot12<- plot(iso2~iso1, data=DS12, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-24), ylim=c(5,11), 
              pch=16, col="black", main="2010")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

plotSO13<- plotSiberObject(my.siber.data13,
                           ax.pad = 2, 
                           hulls = F, community.hulls.args, 
                           ellipses = T, group.ellipses.args,
                           group.hulls = F, group.hull.args,
                           bty = "L",
                           iso.order = c(1,2),
                           xlab = expression({delta}^13*C~'\u2030'),
                           ylab = expression({delta}^15*N~'\u2030'))

c.id<- 13
g.id<- 1

plot13<- plot(iso2~iso1, data=DS13, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-22), ylim=c(5,11), 
              pch=16, col="black", main="2010")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)
#-----
# Merging these

tiff(file="200708.TA.SM.tiff",res=300,width=20,height=20,unit="cm")
par(mfrow=c(1,4))

c.id<- 1
g.id<- 1

plot1<- plot(iso2~iso1, data=DS1, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="black", main="Adelie\nAnvers 2007/08")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

c.id<- 2
g.id<- 1

plot2<- plot(iso2~iso1, data=DS2, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="darkorange", main="Chinstrap\nAnvers 2007/08")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkorange",
                     lty = 1,
                     lwd = 2)

c.id<- 3
g.id<- 1

plot3<- plot(iso2~iso1, data=DS3, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="darkblue", main="Gentoo\nAnvers 2007/08")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkblue",
                     lty = 1,
                     lwd = 2)

c.id<- 10
g.id<- 1

plot10<- plot(iso2~iso1, data=DS10, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
              pch=16, col="darkgray", main="Adelie\nAvian 2007/08")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkgray",
                     lty = 1,
                     lwd = 2)

dev.off()

#-----

tiff(file="200809.TA.SM.tiff",res=300,width=20,height=20,unit="cm")
par(mfrow=c(1,4))

c.id<- 4
g.id<- 1

plot4<- plot(iso2~iso1, data=DS4, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="black", main="Adelie\nAnvers 2008/09")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

c.id<- 5
g.id<- 1

plot5<- plot(iso2~iso1, data=DS5, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="darkorange", main="Chinstrap\nAnvers 2008/09")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkorange",
                     lty = 1,
                     lwd = 2)

c.id<- 6
g.id<- 1

plot6<- plot(iso2~iso1, data=DS6, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="darkblue", main="Gentoo\nAnvers 2008/09")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkblue",
                     lty = 1,
                     lwd = 2)

c.id<- 11
g.id<- 1

plot11<- plot(iso2~iso1, data=DS11, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-22), ylim=c(6,10), 
              pch=16, col="darkgray", main="Adelie\nAvian 2008/09")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkgray",
                     lty = 1,
                     lwd = 2)

dev.off()

#-----

tiff(file="200910.TA.SM.tiff",res=300,width=20,height=20,unit="cm")
par(mfrow=c(1,5))

c.id<- 7
g.id<- 1

plot7<- plot(iso2~iso1, data=DS7, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="black", main="Adelie\nAnvers 2009/10")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "black",
                     lty = 1,
                     lwd = 2)

c.id<- 8
g.id<- 1

plot8<- plot(iso2~iso1, data=DS8, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="darkorange", main="Chinstrap\nAnvers 2009/10")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkorange",
                     lty = 1,
                     lwd = 2)

c.id<- 9
g.id<- 1

plot9<- plot(iso2~iso1, data=DS9, xlab = expression({delta}^13*C~'\u2030'),
             ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
             pch=16, col="darkblue", main="Gentoo\nAnvers 2009/10")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkblue",
                     lty = 1,
                     lwd = 2)

c.id<- 12
g.id<- 1

plot12<- plot(iso2~iso1, data=DS12, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
              pch=16, col="darkgray", main="Adelie\nAvian 2009/10")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "darkgray",
                     lty = 1,
                     lwd = 2)

c.id<- 13
g.id<- 1

plot13<- plot(iso2~iso1, data=DS13, xlab = expression({delta}^13*C~'\u2030'),
              ylab = expression({delta}^15*N~'\u2030'), xlim=c(-29,-21), ylim=c(6,10), 
              pch=16, col="slategray", main="Adelie\nCharcot 2009/10")

coords <- addEllipse(my.siber.data$ML.mu[[c.id]][ , , g.id],
                     my.siber.data$ML.cov[[c.id]][ , , g.id],
                     m = NULL,
                     n = 100,
                     p.interval = 0.95,
                     ci.mean = FALSE,
                     col = "slategray",
                     lty = 1,
                     lwd = 2)

dev.off()



