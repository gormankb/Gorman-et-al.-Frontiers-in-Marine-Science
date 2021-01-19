#-----
# KB Gorman
# Pygoscelis Frontiers Manuscript
# Isotopic mixing model analysis
# 18 Jan, 2020
#-----

# Mixing models using SIAR
# Reset R, clear all objects

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

# SIAR model 1, PAL 2008, NEW 3/19/18, this works with function in code using thinby=15.

data<- read.table('SIAR.PAL.Consumer.5WK.2008.txt', header=TRUE)
source<- read.table('SIAR.PAL.Source1.mndNC.sd.txt', header=TRUE)
corrections<- read.table('SIAR.PAL.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.PAL.Source3.%NC.sd.txt', header=TRUE)

#SIAR.PAL.2008<- siarmcmcdirichletv4(data, source, corrections, concdep) This worked!

#SIAR.PAL.2008<- siarmcmcdirichletv4(data2, source, corrections, concdep, iterations=1000000, burnin=40000, howmany=10000, thinby=10, prior=rep(1,nrow(source)), siardata=list(SHOULDRUN=FALSE)) This is the old model with the code not needed...

SIAR.PAL.2008<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.PAL.2008.plot1<- siarplotdata(SIAR.PAL.2008)

SIAR.PAL.2008.plot2<- siarmatrixplot(SIAR.PAL.2008)
1

SIAR.PAL.2008.plot2<- siarmatrixplot(SIAR.PAL.2008)
2

SIAR.PAL.2008.plot2<- siarmatrixplot(SIAR.PAL.2008)
3

SIAR.PAL.2008.plot3<- siarhistograms(SIAR.PAL.2008)
1
2

SIAR.PAL.2008.plot3<- siarhistograms(SIAR.PAL.2008)
2
2

SIAR.PAL.2008.plot3<- siarhistograms(SIAR.PAL.2008)
3
2

SIAR.PAL.2008.plot4<- siarproportionbygroupplot(SIAR.PAL.2008, prn=TRUE, type="lines")
1
Probability values for Group 1 
95% % lower = 0.67 upper = 0.81 
75% % lower = 0.69 upper = 0.78 
50% % lower = 0.71 upper = 0.76 
Probability values for Group 2 
95% % lower = 0.12 upper = 0.33 
75% % lower = 0.18 upper = 0.3 
50% % lower = 0.21 upper = 0.28 
Probability values for Group 3 
95% % lower = 0 upper = 0.082 
75% % lower = 0.00057 upper = 0.049 
50% % lower = 0.0024 upper = 0.031 

SIAR.PAL.2008.plot4<- siarproportionbygroupplot(SIAR.PAL.2008, prn=TRUE, type="lines")
2
Probability values for Group 1 
95% % lower = 0.7 upper = 0.88 
75% % lower = 0.73 upper = 0.84 
50% % lower = 0.75 upper = 0.82 
Probability values for Group 2 
95% % lower = 0.016 upper = 0.28 
75% % lower = 0.082 upper = 0.25 
50% % lower = 0.13 upper = 0.23 
Probability values for Group 3 
95% % lower = 0 upper = 0.11 
75% % lower = 0.0023 upper = 0.076 
50% % lower = 0.0061 upper = 0.051 

SIAR.PAL.2008.plot4<- siarproportionbygroupplot(SIAR.PAL.2008, prn=TRUE, type="lines")
3
Probability values for Group 1 
95% % lower = 0.52 upper = 0.69 
75% % lower = 0.56 upper = 0.66 
50% % lower = 0.57 upper = 0.63 
Probability values for Group 2 
95% % lower = 0.1 upper = 0.41 
75% % lower = 0.17 upper = 0.35 
50% % lower = 0.21 upper = 0.32 
Probability values for Group 3 
95% % lower = 0.054 upper = 0.21 
75% % lower = 0.085 upper = 0.18 
50% % lower = 0.1 upper = 0.16 


SIAR.PAL.2008.output<- SIAR.PAL.2008$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.PAL.2008.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.PAL.2008.Grp1Sour1<- mean(SIAR.PAL.2008.output[,1])
SIAR.PAL.2008.Grp1Sour1.probval<- "95% % lower = 0.67 upper = 0.81"
SIAR.PAL.2008.Grp1Sour2<- mean(SIAR.PAL.2008.output[,2])
SIAR.PAL.2008.Grp1Sour2.probval<- "95% % lower = 0.12 upper = 0.33"
SIAR.PAL.2008.Grp1Sour3<- mean(SIAR.PAL.2008.output[,3])
SIAR.PAL.2008.Grp1Sour3.probval<- "95% % lower = 0 upper = 0.082"

SIAR.PAL.2008.Grp2Sour1<- mean(SIAR.PAL.2008.output[,6])
SIAR.PAL.2008.Grp2Sour1.probval<- "95% % lower = 0.7 upper = 0.88"
SIAR.PAL.2008.Grp2Sour2<- mean(SIAR.PAL.2008.output[,7])
SIAR.PAL.2008.Grp2Sour2.probval<- "95% % lower = 0.016 upper = 0.28"
SIAR.PAL.2008.Grp2Sour3<- mean(SIAR.PAL.2008.output[,8])
SIAR.PAL.2008.Grp2Sour3.probval<- "95% % lower = 0 upper = 0.11"

SIAR.PAL.2008.Grp3Sour1<- mean(SIAR.PAL.2008.output[,11])
SIAR.PAL.2008.Grp3Sour1.probval<- "95% % lower = 0.52 upper = 0.69"
SIAR.PAL.2008.Grp3Sour2<- mean(SIAR.PAL.2008.output[,12])
SIAR.PAL.2008.Grp3Sour2.probval<- "95% % lower = 0.1 upper = 0.41"
SIAR.PAL.2008.Grp3Sour3<- mean(SIAR.PAL.2008.output[,13])
SIAR.PAL.2008.Grp3Sour3.probval<- "95% % lower = 0.054 upper = 0.21"

SIAR.PAL.2008.mn.AP<- cbind(SIAR.PAL.2008.Grp1Sour1, SIAR.PAL.2008.Grp1Sour1.probval, SIAR.PAL.2008.Grp1Sour2, SIAR.PAL.2008.Grp1Sour2.probval, SIAR.PAL.2008.Grp1Sour3, SIAR.PAL.2008.Grp1Sour3.probval)

SIAR.PAL.2008.mn.CP<- cbind(SIAR.PAL.2008.Grp2Sour1, SIAR.PAL.2008.Grp2Sour1.probval, SIAR.PAL.2008.Grp2Sour2, SIAR.PAL.2008.Grp2Sour2.probval, SIAR.PAL.2008.Grp2Sour3, SIAR.PAL.2008.Grp2Sour3.probval)

SIAR.PAL.2008.mn.GP<- cbind(SIAR.PAL.2008.Grp3Sour1, SIAR.PAL.2008.Grp3Sour1.probval, SIAR.PAL.2008.Grp3Sour2, SIAR.PAL.2008.Grp3Sour2.probval, SIAR.PAL.2008.Grp3Sour3, SIAR.PAL.2008.Grp3Sour3.probval)

SIAR.PAL.2008.mn.all<- as.data.frame(rbind(SIAR.PAL.2008.mn.AP, SIAR.PAL.2008.mn.CP, SIAR.PAL.2008.mn.GP))
names(SIAR.PAL.2008.mn.all)<- c("Esuperb.mn.prop", "95%ProbValue", "Tmacrura.mn.prop", "95%ProbValue", "Eantarctica.mn.prop", "95%ProbValue")
row.names(SIAR.PAL.2008.mn.all)<- c("Adelie", "chinstrap", "gentoo")

write.table(SIAR.PAL.2008.mn.all, file="SIAR.PAL.2008.propmeans&95%probval.csv", col.names=NA, sep=",")

#-----
#PAL 2009

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

# SIAR model 1, PAL 2008, NEW 3/19/18, this works with function in code using thinby=15, but not thinby=10.

data<- read.table('SIAR.PAL.Consumer.5WK.2009.txt', header=TRUE)
source<- read.table('SIAR.PAL.Source1.mndNC.sd.txt', header=TRUE)
corrections<- read.table('SIAR.PAL.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.PAL.Source3.%NC.sd.txt', header=TRUE)

#Load SIAR function

SIAR.PAL.2009<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.PAL.2009.plot1<- siarplotdata(SIAR.PAL.2009)

SIAR.PAL.2009.plot2<- siarmatrixplot(SIAR.PAL.2009)
1

SIAR.PAL.2009.plot2<- siarmatrixplot(SIAR.PAL.2009)
2

SIAR.PAL.2009.plot2<- siarmatrixplot(SIAR.PAL.2009)
3

SIAR.PAL.2009.plot3<- siarhistograms(SIAR.PAL.2009)
1
2

SIAR.PAL.2009.plot3<- siarhistograms(SIAR.PAL.2009)
2
2

SIAR.PAL.2009.plot3<- siarhistograms(SIAR.PAL.2009)
3
2

SIAR.PAL.2009.plot4<- siarproportionbygroupplot(SIAR.PAL.2009, prn=TRUE, type="lines")
1
Probability values for Group 1 
95% % lower = 0.49 upper = 0.62 
75% % lower = 0.52 upper = 0.59 
50% % lower = 0.53 upper = 0.58 
Probability values for Group 2 
95% % lower = 0.3 upper = 0.5 
75% % lower = 0.34 upper = 0.46 
50% % lower = 0.37 upper = 0.44 
Probability values for Group 3 
95% % lower = 0.00098 upper = 0.089 
75% % lower = 0.013 upper = 0.069 
50% % lower = 0.024 upper = 0.059 

SIAR.PAL.2009.plot4<- siarproportionbygroupplot(SIAR.PAL.2009, prn=TRUE, type="lines")
2
Probability values for Group 1 
95% % lower = 0.49 upper = 0.61 
75% % lower = 0.51 upper = 0.59 
50% % lower = 0.53 upper = 0.57 
Probability values for Group 2 
95% % lower = 0.26 upper = 0.47 
75% % lower = 0.3 upper = 0.42 
50% % lower = 0.33 upper = 0.4 
Probability values for Group 3 
95% % lower = 0.032 upper = 0.14 
75% % lower = 0.055 upper = 0.12 
50% % lower = 0.068 upper = 0.1 

SIAR.PAL.2009.plot4<- siarproportionbygroupplot(SIAR.PAL.2009, prn=TRUE, type="lines")
3
Probability values for Group 1 
95% % lower = 0.53 upper = 0.74 
75% % lower = 0.57 upper = 0.69 
50% % lower = 0.6 upper = 0.67 
Probability values for Group 2 
95% % lower = 0.07 upper = 0.35 
75% % lower = 0.13 upper = 0.29 
50% % lower = 0.16 upper = 0.26 
Probability values for Group 3 
95% % lower = 0.09 upper = 0.22 
75% % lower = 0.12 upper = 0.2 
50% % lower = 0.13 upper = 0.18 

SIAR.PAL.2009.output<- SIAR.PAL.2009$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.PAL.2009.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.PAL.2009.Grp1Sour1<- mean(SIAR.PAL.2009.output[,1])
SIAR.PAL.2009.Grp1Sour1.probval<- "95% % lower = 0.49 upper = 0.62"
SIAR.PAL.2009.Grp1Sour2<- mean(SIAR.PAL.2009.output[,2])
SIAR.PAL.2009.Grp1Sour2.probval<- "95% % lower = 0.3 upper = 0.5"
SIAR.PAL.2009.Grp1Sour3<- mean(SIAR.PAL.2009.output[,3])
SIAR.PAL.2009.Grp1Sour3.probval<- "95% % lower = 0.00098 upper = 0.089"

SIAR.PAL.2009.Grp2Sour1<- mean(SIAR.PAL.2009.output[,6])
SIAR.PAL.2009.Grp2Sour1.probval<- "95% % lower = 0.49 upper = 0.61"
SIAR.PAL.2009.Grp2Sour2<- mean(SIAR.PAL.2009.output[,7])
SIAR.PAL.2009.Grp2Sour2.probval<- "95% % lower = 0.26 upper = 0.47"
SIAR.PAL.2009.Grp2Sour3<- mean(SIAR.PAL.2009.output[,8])
SIAR.PAL.2009.Grp2Sour3.probval<- "95% % lower = 0.032 upper = 0.14"

SIAR.PAL.2009.Grp3Sour1<- mean(SIAR.PAL.2009.output[,11])
SIAR.PAL.2009.Grp3Sour1.probval<- "95% % lower = 0.53 upper = 0.7"
SIAR.PAL.2009.Grp3Sour2<- mean(SIAR.PAL.2009.output[,12])
SIAR.PAL.2009.Grp3Sour2.probval<- "95% % lower = 0.07 upper = 0.35"
SIAR.PAL.2009.Grp3Sour3<- mean(SIAR.PAL.2009.output[,13])
SIAR.PAL.2009.Grp3Sour3.probval<- "95% % lower = 0.09 upper = 0.22"

SIAR.PAL.2009.mn.AP<- cbind(SIAR.PAL.2009.Grp1Sour1, SIAR.PAL.2009.Grp1Sour1.probval, SIAR.PAL.2009.Grp1Sour2, SIAR.PAL.2009.Grp1Sour2.probval, SIAR.PAL.2009.Grp1Sour3, SIAR.PAL.2009.Grp1Sour3.probval)

SIAR.PAL.2009.mn.CP<- cbind(SIAR.PAL.2009.Grp2Sour1, SIAR.PAL.2009.Grp2Sour1.probval, SIAR.PAL.2009.Grp2Sour2, SIAR.PAL.2009.Grp2Sour2.probval, SIAR.PAL.2009.Grp2Sour3, SIAR.PAL.2009.Grp2Sour3.probval)

SIAR.PAL.2009.mn.GP<- cbind(SIAR.PAL.2009.Grp3Sour1, SIAR.PAL.2009.Grp3Sour1.probval, SIAR.PAL.2009.Grp3Sour2, SIAR.PAL.2009.Grp3Sour2.probval, SIAR.PAL.2009.Grp3Sour3, SIAR.PAL.2009.Grp3Sour3.probval)

SIAR.PAL.2009.mn.all<- as.data.frame(rbind(SIAR.PAL.2009.mn.AP, SIAR.PAL.2009.mn.CP, SIAR.PAL.2009.mn.GP))
names(SIAR.PAL.2009.mn.all)<- c("Esuperb.mn.prop", "95%ProbValue", "Tmacrura.mn.prop", "95%ProbValue", "Eantarctica.mn.prop", "95%ProbValue")
row.names(SIAR.PAL.2009.mn.all)<- c("Adelie", "chinstrap", "gentoo")

write.table(SIAR.PAL.2009.mn.all, file="SIAR.PAL.2009.propmeans&95%probval.csv", col.names=NA, sep=",")

#-----
#PAL 2010

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

# SIAR model 1, PAL 2008, NEW 3/19/18, this works with function in code using thinby=15, but not thinby=10.

data<- read.table('SIAR.PAL.Consumer.5WK.2010.txt', header=TRUE)
source<- read.table('SIAR.PAL.Source1.mndNC.sd.txt', header=TRUE)
corrections<- read.table('SIAR.PAL.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.PAL.Source3.%NC.sd.txt', header=TRUE)

#Load SIAR function

SIAR.PAL.2010<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.PAL.2010.plot1<- siarplotdata(SIAR.PAL.2010)

SIAR.PAL.2010.plot2<- siarmatrixplot(SIAR.PAL.2010)
1

SIAR.PAL.2010.plot2<- siarmatrixplot(SIAR.PAL.2010)
2

SIAR.PAL.2010.plot2<- siarmatrixplot(SIAR.PAL.2010)
3

SIAR.PAL.2010.plot3<- siarhistograms(SIAR.PAL.2010)
1
2

SIAR.PAL.2010.plot3<- siarhistograms(SIAR.PAL.2010)
2
2

SIAR.PAL.2010.plot3<- siarhistograms(SIAR.PAL.2010)
3
2

SIAR.PAL.2010.plot4<- siarproportionbygroupplot(SIAR.PAL.2010, prn=TRUE, type="lines")
1
Probability values for Group 1 
95% % lower = 0.59 upper = 0.76 
75% % lower = 0.62 upper = 0.72 
50% % lower = 0.64 upper = 0.7 
Probability values for Group 2 
95% % lower = 0.021 upper = 0.3 
75% % lower = 0.079 upper = 0.25 
50% % lower = 0.12 upper = 0.22 
Probability values for Group 3 
95% % lower = 0.086 upper = 0.24 
75% % lower = 0.12 upper = 0.21 
50% % lower = 0.13 upper = 0.19 

SIAR.PAL.2010.plot4<- siarproportionbygroupplot(SIAR.PAL.2010, prn=TRUE, type="lines")
2
Probability values for Group 1 
95% % lower = 0.54 upper = 0.71 
75% % lower = 0.57 upper = 0.67 
50% % lower = 0.59 upper = 0.65 
Probability values for Group 2 
95% % lower = 0.073 upper = 0.37 
75% % lower = 0.14 upper = 0.31 
50% % lower = 0.18 upper = 0.27 
Probability values for Group 3 
95% % lower = 0.076 upper = 0.23 
75% % lower = 0.11 upper = 0.2 
50% % lower = 0.13 upper = 0.18 

SIAR.PAL.2010.plot4<- siarproportionbygroupplot(SIAR.PAL.2010, prn=TRUE, type="lines")
3
Probability values for Group 1 
95% % lower = 0.64 upper = 0.76 
75% % lower = 0.67 upper = 0.74 
50% % lower = 0.68 upper = 0.73 
Probability values for Group 2 
95% % lower = 0 upper = 0.14 
75% % lower = 0.00035 upper = 0.082 
50% % lower = 0.0033 upper = 0.051 
Probability values for Group 3 
95% % lower = 0.18 upper = 0.3 
75% % lower = 0.21 upper = 0.28 
50% % lower = 0.23 upper = 0.27 

SIAR.PAL.2010.output<- SIAR.PAL.2010$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.PAL.2010.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.PAL.2010.Grp1Sour1<- mean(SIAR.PAL.2010.output[,1])
SIAR.PAL.2010.Grp1Sour1.probval<- "95% % lower = 0.59 upper = 0.76"
SIAR.PAL.2010.Grp1Sour2<- mean(SIAR.PAL.2010.output[,2])
SIAR.PAL.2010.Grp1Sour2.probval<- "95% % lower = 0.021 upper = 0.3"
SIAR.PAL.2010.Grp1Sour3<- mean(SIAR.PAL.2010.output[,3])
SIAR.PAL.2010.Grp1Sour3.probval<- "95% % lower = 0.086 upper = 0.24"

SIAR.PAL.2010.Grp2Sour1<- mean(SIAR.PAL.2010.output[,6])
SIAR.PAL.2010.Grp2Sour1.probval<- "95% % lower = 0.54 upper = 0.71"
SIAR.PAL.2010.Grp2Sour2<- mean(SIAR.PAL.2010.output[,7])
SIAR.PAL.2010.Grp2Sour2.probval<- "95% % lower = 0.073 upper = 0.37"
SIAR.PAL.2010.Grp2Sour3<- mean(SIAR.PAL.2010.output[,8])
SIAR.PAL.2010.Grp2Sour3.probval<- "95% % lower = 0.076 upper = 0.23"

SIAR.PAL.2010.Grp3Sour1<- mean(SIAR.PAL.2010.output[,11])
SIAR.PAL.2010.Grp3Sour1.probval<- "95% % lower = 0.64 upper = 0.76"
SIAR.PAL.2010.Grp3Sour2<- mean(SIAR.PAL.2010.output[,12])
SIAR.PAL.2010.Grp3Sour2.probval<- "95% % lower = 0 upper = 0.14"
SIAR.PAL.2010.Grp3Sour3<- mean(SIAR.PAL.2010.output[,13])
SIAR.PAL.2010.Grp3Sour3.probval<- "95% % lower = 0.18 upper = 0.3"

SIAR.PAL.2010.mn.AP<- cbind(SIAR.PAL.2010.Grp1Sour1, SIAR.PAL.2010.Grp1Sour1.probval, SIAR.PAL.2010.Grp1Sour2, SIAR.PAL.2010.Grp1Sour2.probval, SIAR.PAL.2010.Grp1Sour3, SIAR.PAL.2010.Grp1Sour3.probval)

SIAR.PAL.2010.mn.CP<- cbind(SIAR.PAL.2010.Grp2Sour1, SIAR.PAL.2010.Grp2Sour1.probval, SIAR.PAL.2010.Grp2Sour2, SIAR.PAL.2010.Grp2Sour2.probval, SIAR.PAL.2010.Grp2Sour3, SIAR.PAL.2010.Grp2Sour3.probval)

SIAR.PAL.2010.mn.GP<- cbind(SIAR.PAL.2010.Grp3Sour1, SIAR.PAL.2010.Grp3Sour1.probval, SIAR.PAL.2010.Grp3Sour2, SIAR.PAL.2010.Grp3Sour2.probval, SIAR.PAL.2010.Grp3Sour3, SIAR.PAL.2010.Grp3Sour3.probval)

SIAR.PAL.2010.mn.all<- as.data.frame(rbind(SIAR.PAL.2010.mn.AP, SIAR.PAL.2010.mn.CP, SIAR.PAL.2010.mn.GP))
names(SIAR.PAL.2010.mn.all)<- c("Esuperb.mn.prop", "95%ProbValue", "Tmacrura.mn.prop", "95%ProbValue", "Eantarctica.mn.prop", "95%ProbValue")
row.names(SIAR.PAL.2010.mn.all)<- c("Adelie", "chinstrap", "gentoo")

write.table(SIAR.PAL.2010.mn.all, file="SIAR.PAL.2010.propmeans&95%probval.csv", col.names=NA, sep=",")

#-----
# Avian, 2008

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

data<- read.table('SIAR.AVI.Consumer.5WK.2008.txt', header=TRUE)
source<- read.table('SIAR.AVI.Source1.mndNC.sd.txt', header=TRUE)
corrections<- read.table('SIAR.AVI.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.AVI.Source3.%NC.sd.txt', header=TRUE)

#Load SIAR function

SIAR.AVI.2008<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.AVI.2008.plot1<- siarplotdata(SIAR.AVI.2008)

SIAR.AVI.2008.plot2<- siarmatrixplot(SIAR.AVI.2008)
1


SIAR.AVI.2008.plot3<- siarhistograms(SIAR.AVI.2008)
1
2


SIAR.AVI.2008.plot4<- siarproportionbygroupplot(SIAR.AVI.2008, prn=TRUE, type="lines")
1

Probability values for Group 1 
95% % lower = 0.085 upper = 0.23 
75% % lower = 0.12 upper = 0.21 
50% % lower = 0.14 upper = 0.19 
Probability values for Group 2 
95% % lower = 0.48 upper = 0.65 
75% % lower = 0.51 upper = 0.61 
50% % lower = 0.53 upper = 0.59 
Probability values for Group 3 
95% % lower = 0.056 upper = 0.29 
75% % lower = 0.11 upper = 0.25 
50% % lower = 0.14 upper = 0.22 
Probability values for Group 4 
95% % lower = 0 upper = 0.14 
75% % lower = 0.0024 upper = 0.092 
50% % lower = 0.0071 upper = 0.062 
Probability values for Group 5 
95% % lower = 0 upper = 0.09 
75% % lower = 0.00021 upper = 0.051 
50% % lower = 0.0019 upper = 0.031 


SIAR.AVI.2008.output<- SIAR.AVI.2008$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.AVI.2008.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.AVI.2008.Grp1Sour1<- mean(SIAR.AVI.2008.output[,1])
SIAR.AVI.2008.Grp1Sour1.probval<- "95% % lower = 0.085 upper = 0.23"
SIAR.AVI.2008.Grp1Sour2<- mean(SIAR.AVI.2008.output[,2])
SIAR.AVI.2008.Grp1Sour2.probval<- "95% % lower = 0.48 upper = 0.65"
SIAR.AVI.2008.Grp1Sour3<- mean(SIAR.AVI.2008.output[,3])
SIAR.AVI.2008.Grp1Sour3.probval<- "95% % lower = 0.056 upper = 0.29"
SIAR.AVI.2008.Grp1Sour4<- mean(SIAR.AVI.2008.output[,4])
SIAR.AVI.2008.Grp1Sour4.probval<- "95% % lower = 0 upper = 0.14"
SIAR.AVI.2008.Grp1Sour5<- mean(SIAR.AVI.2008.output[,5])
SIAR.AVI.2008.Grp1Sour5.probval<- "95% % lower = 0 upper = 0.09"

SIAR.AVI.2008.mn.all<- as.data.frame(cbind(SIAR.AVI.2008.Grp1Sour1, SIAR.AVI.2008.Grp1Sour1.probval, SIAR.AVI.2008.Grp1Sour2, SIAR.AVI.2008.Grp1Sour2.probval, SIAR.AVI.2008.Grp1Sour3, SIAR.AVI.2008.Grp1Sour3.probval, SIAR.AVI.2008.Grp1Sour4, SIAR.AVI.2008.Grp1Sour4.probval, SIAR.AVI.2008.Grp1Sour5, SIAR.AVI.2008.Grp1Sour5.probval))
names(SIAR.AVI.2008.mn.all)<- c("Ecrystal.mn.prop", "95%ProbValue", "Esuperb.mn.prop", "95%ProbValue", "Tmacrura.mn.prop", "95%ProbValue", "Eantarctica.mn.prop", "95%ProbValue", "Pantarcticum.mn.prob", "95%ProbValue")
row.names(SIAR.AVI.2008.mn.all)<- c("Adelie")

write.table(SIAR.AVI.2008.mn.all, file="SIAR.AVI.2008.propmeans&95%probval.csv", col.names=NA, sep=",")

#-----
# Avian, 2009

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

data<- read.table('SIAR.AVI.Consumer.5WK.2009.txt', header=TRUE)
source<- read.table('SIAR.AVI.Source1.mndNC.sd.txt', header=TRUE)
corrections<- read.table('SIAR.AVI.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.AVI.Source3.%NC.sd.txt', header=TRUE)

#Load SIAR function

SIAR.AVI.2009<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.AVI.2009.plot1<- siarplotdata(SIAR.AVI.2009)

SIAR.AVI.2009.plot2<- siarmatrixplot(SIAR.AVI.2009)
1



SIAR.AVI.2009.plot3<- siarhistograms(SIAR.AVI.2009)
1
2


SIAR.AVI.2009.plot4<- siarproportionbygroupplot(SIAR.AVI.2009, prn=TRUE, type="lines")
1
Probability values for Group 1 
95% % lower = 0.22 upper = 0.38 
75% % lower = 0.26 upper = 0.35 
50% % lower = 0.28 upper = 0.34 
Probability values for Group 2 
95% % lower = 0.56 upper = 0.67 
75% % lower = 0.58 upper = 0.65 
50% % lower = 0.6 upper = 0.63 
Probability values for Group 3 
95% % lower = 0 upper = 0.081 
75% % lower = 0.00091 upper = 0.05 
50% % lower = 0.0028 upper = 0.032 
Probability values for Group 4 
95% % lower = 0 upper = 0.071 
75% % lower = 0 upper = 0.039 
50% % lower = 0.0013 upper = 0.023 
Probability values for Group 5 
95% % lower = 0 upper = 0.054 
75% % lower = 0 upper = 0.028 
50% % lower = 0.0008 upper = 0.016 

SIAR.AVI.2009.output<- SIAR.AVI.2009$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.AVI.2009.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.AVI.2009.Grp1Sour1<- mean(SIAR.AVI.2009.output[,1])
SIAR.AVI.2009.Grp1Sour1.probval<- "95% % lower = 0.22 upper = 0.38"
SIAR.AVI.2009.Grp1Sour2<- mean(SIAR.AVI.2009.output[,2])
SIAR.AVI.2009.Grp1Sour2.probval<- "95% % lower = 0.56 upper = 0.67"
SIAR.AVI.2009.Grp1Sour3<- mean(SIAR.AVI.2009.output[,3])
SIAR.AVI.2009.Grp1Sour3.probval<- "95% % lower = 0 upper = 0.081"
SIAR.AVI.2009.Grp1Sour4<- mean(SIAR.AVI.2009.output[,4])
SIAR.AVI.2009.Grp1Sour4.probval<- "95% % lower = 0 upper = 0.071"
SIAR.AVI.2009.Grp1Sour5<- mean(SIAR.AVI.2009.output[,5])
SIAR.AVI.2009.Grp1Sour5.probval<- "95% % lower = 0 upper = 0.054"

SIAR.AVI.2009.mn.all<- as.data.frame(cbind(SIAR.AVI.2009.Grp1Sour1, SIAR.AVI.2009.Grp1Sour1.probval, SIAR.AVI.2009.Grp1Sour2, SIAR.AVI.2009.Grp1Sour2.probval, SIAR.AVI.2009.Grp1Sour3, SIAR.AVI.2009.Grp1Sour3.probval, SIAR.AVI.2009.Grp1Sour4, SIAR.AVI.2009.Grp1Sour4.probval, SIAR.AVI.2009.Grp1Sour5, SIAR.AVI.2009.Grp1Sour5.probval))
names(SIAR.AVI.2009.mn.all)<- c("Ecrystal.mn.prop", "95%ProbValue", "Esuperb.mn.prop", "95%ProbValue", "Tmacrura.mn.prop", "95%ProbValue", "Eantarctica.mn.prop", "95%ProbValue", "Pantarcticum.mn.prob", "95%ProbValue")
row.names(SIAR.AVI.2009.mn.all)<- c("Adelie")

write.table(SIAR.AVI.2009.mn.all, file="SIAR.AVI.2009.propmeans&95%probval.csv", col.names=NA, sep=",")

#-----
# Avian, 2010

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

data<- read.table('SIAR.AVI.Consumer.5WK.2010.txt', header=TRUE)
source<- read.table('SIAR.AVI.Source1.mndNC.sd-2.txt', header=TRUE)
corrections<- read.table('SIAR.AVI.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.AVI.Source3.%NC.sd-2.txt', header=TRUE)

#Load SIAR function

SIAR.AVI.2010<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.AVI.2010.plot1<- siarplotdata(SIAR.AVI.2010)

SIAR.AVI.2010.plot2<- siarmatrixplot(SIAR.AVI.2010)
1



SIAR.AVI.2010.plot3<- siarhistograms(SIAR.AVI.2010)
1
2


SIAR.AVI.2010.plot4<- siarproportionbygroupplot(SIAR.AVI.2010, prn=TRUE, type="lines")
1
Probability values for Group 1 
95% % lower = 0.015 upper = 0.14 
75% % lower = 0.042 upper = 0.11 
50% % lower = 0.057 upper = 0.099 
Probability values for Group 2 
95% % lower = 0.52 upper = 0.7 
75% % lower = 0.55 upper = 0.65 
50% % lower = 0.57 upper = 0.63 
Probability values for Group 3 
95% % lower = 0.14 upper = 0.4 
75% % lower = 0.2 upper = 0.36 
50% % lower = 0.24 upper = 0.33 
Probability values for Group 4 
95% % lower = 0 upper = 0.071 
75% % lower = 0 upper = 0.037 
50% % lower = 0.0011 upper = 0.022 
Probability values for Group 5 
95% % lower = 0 upper = 0.046 
75% % lower = 0 upper = 0.023 
50% % lower = 0.00072 upper = 0.014 

# Run 2, 4/7/18
Probability values for Group 1 
95% % lower = 0.041 upper = 0.17 
75% % lower = 0.071 upper = 0.15 
50% % lower = 0.088 upper = 0.13 
Probability values for Group 2 
95% % lower = 0.55 upper = 0.74 
75% % lower = 0.58 upper = 0.7 
50% % lower = 0.6 upper = 0.67 
Probability values for Group 3 
95% % lower = 0.042 upper = 0.32 
75% % lower = 0.11 upper = 0.28 
50% % lower = 0.15 upper = 0.25 
Probability values for Group 4 
95% % lower = 0 upper = 0.092 
75% % lower = 0.000025 upper = 0.052 
50% % lower = 0.0018 upper = 0.032 
Probability values for Group 5 
95% % lower = 0 upper = 0.071 
75% % lower = 0 upper = 0.039 
50% % lower = 0.0012 upper = 0.023 


SIAR.AVI.2010.output<- SIAR.AVI.2010$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.AVI.2010.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.AVI.2010.Grp1Sour1<- mean(SIAR.AVI.2010.output[,1])
SIAR.AVI.2010.Grp1Sour1.probval<- "95% % lower = 0.041 upper = 0.17"
SIAR.AVI.2010.Grp1Sour2<- mean(SIAR.AVI.2010.output[,2])
SIAR.AVI.2010.Grp1Sour2.probval<- "95% % lower = 0.55 upper = 0.74"
SIAR.AVI.2010.Grp1Sour3<- mean(SIAR.AVI.2010.output[,3])
SIAR.AVI.2010.Grp1Sour3.probval<- "95% % lower = 0.042 upper = 0.32"
SIAR.AVI.2010.Grp1Sour4<- mean(SIAR.AVI.2010.output[,4])
SIAR.AVI.2010.Grp1Sour4.probval<- "95% % lower = 0 upper = 0.092"
SIAR.AVI.2010.Grp1Sour5<- mean(SIAR.AVI.2010.output[,5])
SIAR.AVI.2010.Grp1Sour5.probval<- "95% % lower = 0 upper = 0.071"

SIAR.AVI.2010.mn.all<- as.data.frame(cbind(SIAR.AVI.2010.Grp1Sour1, SIAR.AVI.2010.Grp1Sour1.probval, SIAR.AVI.2010.Grp1Sour2, SIAR.AVI.2010.Grp1Sour2.probval, SIAR.AVI.2010.Grp1Sour3, SIAR.AVI.2010.Grp1Sour3.probval, SIAR.AVI.2010.Grp1Sour4, SIAR.AVI.2010.Grp1Sour4.probval, SIAR.AVI.2010.Grp1Sour5, SIAR.AVI.2010.Grp1Sour5.probval))
names(SIAR.AVI.2010.mn.all)<- c("Ecrystal.mn.prop", "95%ProbValue", "Esuperb.mn.prop", "95%ProbValue", "Tmacrura.mn.prop", "95%ProbValue", "Eantarctica.mn.prop", "95%ProbValue", "Pantarcticum.mn.prob", "95%ProbValue")
row.names(SIAR.AVI.2010.mn.all)<- c("Adelie")

write.table(SIAR.AVI.2010.mn.all, file="SIAR.AVI.2010.propmeans&95%probval.csv", col.names=NA, sep=",")

#-----
# CHA, -2

rm(list=ls())

library(devtools)
install_github("andrewljackson/siar@master", build_vingettes = TRUE)
library(siar)

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

data<- read.table('SIAR.CHA.Consumer.5WK.2010.txt', header=TRUE)
source<- read.table('SIAR.CHA.Source1.mndNC.sd-2.txt', header=TRUE)
corrections<- read.table('SIAR.CHA.Source2.DFdNC.sd.txt', header=TRUE)
concdep<- read.table('SIAR.CHA.Source3.%NC.sd-2.txt', header=TRUE)

#Load SIAR function

SIAR.CHA.2010<- siarmcmcdirichletv4(data, source, corrections, concdep)

# (2000000-50000)/15=130000 posterior draws

SIAR.CHA.2010.plot1<- siarplotdata(SIAR.CHA.2010)

SIAR.CHA.2010.plot2<- siarmatrixplot(SIAR.CHA.2010)
1



SIAR.CHA.2010.plot3<- siarhistograms(SIAR.CHA.2010)
1
2


SIAR.CHA.2010.plot4<- siarproportionbygroupplot(SIAR.CHA.2010, prn=TRUE, type="lines")
1
Probability values for Group 1 
95% % lower = 0.39 upper = 0.67 
75% % lower = 0.45 upper = 0.62 
50% % lower = 0.5 upper = 0.59 
Probability values for Group 2 
95% % lower = 0.22 upper = 0.45 
75% % lower = 0.26 upper = 0.4 
50% % lower = 0.29 upper = 0.37 
Probability values for Group 3 
95% % lower = 0.023 upper = 0.24 
75% % lower = 0.071 upper = 0.2 
50% % lower = 0.099 upper = 0.18 
#run 2
Probability values for Group 1 
95% % lower = 0.37 upper = 0.67 
75% % lower = 0.44 upper = 0.62 
50% % lower = 0.48 upper = 0.59 
Probability values for Group 2 
95% % lower = 0.23 upper = 0.45 
75% % lower = 0.27 upper = 0.4 
50% % lower = 0.29 upper = 0.37 
Probability values for Group 3 
95% % lower = 0.027 upper = 0.25 
75% % lower = 0.077 upper = 0.21 
50% % lower = 0.1 upper = 0.18 

SIAR.CHA.2010.output<- SIAR.CHA.2010$output # Don't want to call this, output too big. Just examine with fix fnx.

fix(SIAR.CHA.2010.output) # Couldnt get this to save file!


# Means of groups proportions using output
SIAR.CHA.2010.Grp1Sour1<- mean(SIAR.CHA.2010.output[,1])
SIAR.CHA.2010.Grp1Sour1.probval<- "95% % lower = 0.37 upper = 0.67"
SIAR.CHA.2010.Grp1Sour2<- mean(SIAR.CHA.2010.output[,2])
SIAR.CHA.2010.Grp1Sour2.probval<- "95% % lower = 0.23 upper = 0.4"
SIAR.CHA.2010.Grp1Sour3<- mean(SIAR.CHA.2010.output[,3])
SIAR.CHA.2010.Grp1Sour3.probval<- "95% % lower = 0.027 upper = 0.25"

SIAR.CHA.2010.mn.all<- as.data.frame(cbind(SIAR.CHA.2010.Grp1Sour1, SIAR.CHA.2010.Grp1Sour1.probval, SIAR.CHA.2010.Grp1Sour2, SIAR.CHA.2010.Grp1Sour2.probval, SIAR.CHA.2010.Grp1Sour3, SIAR.CHA.2010.Grp1Sour3.probval))
names(SIAR.CHA.2010.mn.all)<- c("Ecrystal.mn.prop", "95%ProbValue", "Esuperb.mn.prop", "95%ProbValue", "Pantarcticum.mn.prob", "95%ProbValue")
row.names(SIAR.CHA.2010.mn.all)<- c("Adelie")

write.table(SIAR.CHA.2010.mn.all, file="SIAR.CHA.2010.propmeans&95%probval.csv", col.names=NA, sep=",")


#-----
# Plotting Mixing Models

rm(list=ls())
library(ggplot2)
library(gridExtra)

#-----
# Load Data

#-----
setwd("/Users/kgorman/Documents/KB Gorman Files/R Data")
getwd()
list.files()

DS.PAL.2008<- read.csv("SIAR.PAL.2008.ggplot.mnprop95prob.csv")
DS.PAL.2009<- read.csv("SIAR.PAL.2009.ggplot.mnprop95prob.csv")
DS.PAL.2010<- read.csv("SIAR.PAL.2010.ggplot.mnprop95prob.csv")

tiff(file="PAL.SIAR18.tiff",res=150,width=20,height=20,unit="cm")

plot1<- ggplot(data=DS.PAL.2008, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2007/08", y="Proportion of Diet (%)") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("salmon", "orangered", "dodgerblue4")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

plot2<- ggplot(data=DS.PAL.2009, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2008/09", y="") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("salmon", "orangered", "dodgerblue4")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

plot3<- ggplot(data=DS.PAL.2010, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2009/10", y="") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("salmon", "orangered", "dodgerblue4")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

multi.plot<- grid.arrange(plot1, plot2, plot3, ncol=3)
multi.plot

dev.off()

#-----
DS.AVI.2008<- read.csv("SIAR.AVI.2008.ggplot.mnprop95prob.csv")
DS.AVI.2009<- read.csv("SIAR.AVI.2009.ggplot.mnprop95prob.csv")
DS.AVI.2010<- read.csv("SIAR.AVI.2010.ggplot.mnprop95prob.csv")

tiff(file="AVI.SIAR18-2.tiff",res=150,width=20,height=20,unit="cm")

plot1<- ggplot(data=DS.AVI.2008, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2007/08", y="Proportion of Diet (%)") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("pink", "salmon", "orangered", "dodgerblue4", "slategrey")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

plot2<- ggplot(data=DS.AVI.2009, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2008/09", y="") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("pink", "salmon", "orangered", "dodgerblue4", "slategrey")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

plot3<- ggplot(data=DS.AVI.2010, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2009/10", y="") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("pink", "salmon", "orangered", "dodgerblue4", "slategrey")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

multi.plot<- grid.arrange(plot1, plot2, plot3, ncol=3)
multi.plot

dev.off()

#-----

DS.CHA.2010<- read.csv("SIAR.CHA.2010.ggplot.mnprop95prob.csv")

tiff(file="CHA.SIAR18.tiff",res=150,width=20,height=20,unit="cm")

plot1<- ggplot(data=DS.CHA.2010, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2007/08", y="") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("pink", "salmon", "slategrey")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

plot2<- ggplot(data=DS.CHA.2010, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2008/09", y="") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("pink", "salmon", "slategrey")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

plot3<- ggplot(data=DS.CHA.2010, aes(fill=Prey, y=Proportion, x=Species)) + geom_bar(position="dodge", stat="identity") + labs(x="2009/10", y="Proportion of Diet (%)") + ylim(0,1) + theme(axis.title.x=element_text(vjust=-.5, size=16)) + theme(axis.title.y=element_text(angle=90, vjust=2, size=18)) + theme(axis.text.x=element_text(size=10)) + theme(axis.text.y=element_text(size=10)) + scale_fill_manual(values=c("pink", "salmon", "slategrey")) + theme(plot.margin=unit(c(.5,.25,.5,.5), "cm")) + theme(legend.position="none") + geom_errorbar(aes(ymin = Proportion-X95.lower.2, ymax = Proportion+X95.upper.2), width=0.2, position=position_dodge(0.9))

multi.plot<- grid.arrange(plot1, plot2, plot3, ncol=3)
multi.plot

dev.off()