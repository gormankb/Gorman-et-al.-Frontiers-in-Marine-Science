#-----
# KB Gorman
# Pygoscelis Frontiers Manuscript
# Isotopic Dietary Analysis (IDA)
# 30 Nov, 2020
#-----
# Variation in C/N between species analyses. Adults, IDA analysis 1.

# Reset R, clear all objects.
rm(list=ls())

# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see what files are available.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the data to be used.
DataFileName<- "Ads.INS.DataSet.Sub4.csv" #Ads C
DataFileName<- "Ads.INS.DataSet.Sub3.csv" #Ads N

# Identify the models to be run.
ModelFileName<-"Ads.IsoNicheSpace.CandSet.C.csv"
ModelFileName<-"Ads.IsoNicheSpace.CandSet.N.csv"

#-----
# Import the data. Be sure to check the structure of the data and that R is reading continuous and categorical variables appropriately.
DataSet<-read.csv(DataFileName)

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

# Change parameter structures if needed.
DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

# Import candidate model set and convert to a character vector.
ModelSet.df<-read.csv(ModelFileName)  
ModelSet<-as.character(ModelSet.df$Model)

# View candidate models.
ModelSet

#-----
# Calculate all combinations of predictor variables to check candidate model set. Note, any candidate set may not include all combos of variables, but this will allow for checking that all combos have at least been considered (i.e., didnt forget some combo).

globalmodel<- lm(delta.15.N~Spp*Year.2*Spp:Year.2, data=DataSet)

AllCombCoefName<-names(globalmodel$coef)

# View all combinations of predictor variables.
AllCombCoefName

#-----
# Specify names of main explanatory parameters and associated SEs in the ModelSet (including interactions). Need to use : for interactions. !!Make sure that when creating any subsequent versions of CoefName in the code below that it follows the same order as is listed here. Further, when a candidate set has interaction models where (a) both main parameters are not included in the model, or where (b) the main effect is a categorical variable, you must include columns labeled here for all appropriate output!!
CoefName<- c("Intercept", "Intercept.SE", "SppCHPE", "SppCHPE.SE", "SppGEPE", "SppGEPE.SE", "Year.22", "Year.22.SE", "Year.23", "Year.23.SE", "SppCHPE:Year.22", "SppCHPE:Year.22.SE", "SppCHPE:Year.23", "SppCHPE:Year.23.SE", "SppGEPE:Year.22", "SppGEPE:Year.22.SE", "SppGEPE:Year.23", "SppGEPE:Year.23.SE")

# View predictor variables and associated SEs.
CoefName

#-----
# AIC function produces calculations for AIC Model Matrix that is used in conjunction with the ModelSet.
calculate.AIC<- function(AIC.Table, ModelSet) {
  
  deltaAIC.c<- AIC.Table$AIC.c-min(AIC.Table$AIC.c)
  lik.dAIC.c<- exp(-deltaAIC.c/2)
  AIC.c.W<- lik.dAIC.c/(sum(lik.dAIC.c))
  AICFinalMatrix<- data.frame(AICModelMatrix, deltaAIC.c, lik.dAIC.c, AIC.c.W)
}

#-----
# Create matrix to hold output from models.

# I. AIC Model Matrix, ncol=29 includes models, all coef names listed above, and column headings specified in the matrix below.
AICModelMatrix<- as.data.frame(matrix(NA, nrow=length(ModelSet), ncol=length(CoefName)+11, dimnames=list(c(1:length(ModelSet)), c("Models", CoefName, "N.Obs", "k", "EDF", "RMSE", "SSE", "logLik", "-2logLik", "mul.r.squared", "AIC", "AIC.c"))))

# View AIC Model Matrix.
head(AICModelMatrix)

#-----
# Loop for calculating model output and filling AIC Model Matrix. According to Burnham and Anderson, for least squares model fitting, K = total number of estimated regression parameters, including the intercept, and residual variation. In the case of R output, coef(model) includes an estimate of the intercept, thus, length(m$coef)+1 = K.

ModelOutput<- list()

for(i in 1:length(ModelSet)){
  
  ModelOutput[[i]]<- lm(as.formula(ModelSet[i]), data=DataSet)
  m<- ModelOutput[[i]]
  N.Obs<- nrow(DataSet)
  k<- length(m$coef)+1
  EDF<- N.Obs-k
  AIC<- AIC(m)
  AIC.c<- AIC+(2*k*(k+1))/(N.Obs-k-1)
  AICModelMatrix[i,"Models"]<- ModelSet[i]
  AICModelMatrix[i,"Intercept"]<- coef(m)["(Intercept)"]
  AICModelMatrix[i,"Intercept.SE"]<- summary(m)$coef["(Intercept)",2]
  AICModelMatrix[i,"SppCHPE"]<- coef(m)["SppCHPE"]
  AICModelMatrix[i,"SppCHPE.SE"]<- summary(m)$coef[,2]["SppCHPE"]
  AICModelMatrix[i,"SppGEPE"]<- coef(m)["SppGEPE"]
  AICModelMatrix[i,"SppGEPE.SE"]<- summary(m)$coef[,2]["SppGEPE"]
  AICModelMatrix[i,"Year.22"]<- coef(m)["Year.22"]
  AICModelMatrix[i,"Year.22.SE"]<- summary(m)$coef[,2]["Year.22"]
  AICModelMatrix[i,"Year.23"]<- coef(m)["Year.23"]
  AICModelMatrix[i,"Year.23.SE"]<- summary(m)$coef[,2]["Year.23"]
  AICModelMatrix[i,"SppCHPE:Year.22"]<- coef(m)["SppCHPE:Year.22"]
  AICModelMatrix[i,"SppCHPE:Year.22.SE"]<- summary(m)$coef[,2]["SppCHPE:Year.22"]
  AICModelMatrix[i,"SppCHPE:Year.23"]<- coef(m)["SppCHPE:Year.23"]
  AICModelMatrix[i,"SppCHPE:Year.23.SE"]<- summary(m)$coef[,2]["SppCHPE:Year.23"]
  AICModelMatrix[i,"SppGEPE:Year.22"]<- coef(m)["SppGEPE:Year.22"]
  AICModelMatrix[i,"SppGEPE:Year.22.SE"]<- summary(m)$coef[,2]["SppGEPE:Year.22"]
  AICModelMatrix[i,"SppGEPE:Year.23"]<- coef(m)["SppGEPE:Year.23"]
  AICModelMatrix[i,"SppGEPE:Year.23.SE"]<- summary(m)$coef[,2]["SppGEPE:Year.23"]
  AICModelMatrix[i,"N.Obs"]<- N.Obs
  AICModelMatrix[i,"k"]<- k 
  AICModelMatrix[i,"EDF"]<- EDF
  AICModelMatrix[i,"RMSE"]<- summary(m)$sigma
  AICModelMatrix[i,"SSE"]<- anova(m)["Residuals","Sum Sq"]
  AICModelMatrix[i,"logLik"]<- logLik(m)
  AICModelMatrix[i,"-2logLik"]<- -2*logLik(m)
  AICModelMatrix[i,"mul.r.squared"]<- summary(m)$r.squared
  AICModelMatrix[i,"AIC"]<- AIC
  AICModelMatrix[i,"AIC.c"]<- AIC.c
}

#-----
# Calculate deltaAICc, likdAICc, and AIC weights.
AIC.Output<-calculate.AIC(AICModelMatrix,as.character(ModelSet))

#-----
# View AIC Model Matrix.
print(AIC.Output)

# Write AIC Final Matrix output to .csv files.
write.table(AIC.Output, file="AICFinalMatrix.csv", col.names=NA, sep=",")

# Write model output to .csv files.
sink(paste("Model Summaries_",out="",".doc",sep=""))

for (i in 1:length(ModelSet)) {
  
  print(paste("MODEL ", i, sep = ""))
  print(ModelSet[[i]])
  print(summary(ModelOutput[[i]]))
  print(paste("-------------------------------------------------------------------"))
}

sink()

#-----
# Begin weighted and summed AIC calcs. Save previous file with AICFinalMatrix output, 
# excluding column A (Model #) and N.Obs through the lik.dAIC.c columns, 
# as weighted AIC Model Matrix Worksheet (W.AICFinalMatrixWS.csv).

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see whats available. Make sure that new file created above is listed.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the new data WS to be used.
DataFileName2<- "W.AICFinalMatrixWS.csv"

# Import data.
DataSet2<-read.csv(DataFileName2)

# View DataSet. Note that all columns, with the exception of Models, should be listed as num.
head(DataSet2)
str(DataSet2)

#-----
# Specify a vector for DataSet2 Parameters. These will be used to define the parameters that will used to calculate Parameter likelihoods and W.ParaEsts that are multiplied by the model W. !!For interactions, be sure to specify them with . and not : as this is out they will be uploaded in DataSet2!!

# First recall CoefName and just exclude the .SE names in CoefName2 below.
CoefName<- c("Intercept", "Intercept.SE", "SppCHPE", "SppCHPE.SE", "SppGEPE", "SppGEPE.SE", "Year.22", "Year.22.SE", "Year.23", "Year.23.SE", "SppCHPE:Year.22", "SppCHPE:Year.22.SE", "SppCHPE:Year.23", "SppCHPE:Year.23.SE", "SppGEPE:Year.22", "SppGEPE:Year.22.SE", "SppGEPE:Year.23", "SppGEPE:Year.23.SE")

# View CoefName.
CoefName

CoefName2<- c("Intercept", "SppCHPE", "SppGEPE", "Year.22", "Year.23", "SppCHPE.Year.22", "SppCHPE.Year.23", "SppGEPE.Year.22", "SppGEPE.Year.23")

# View CoefName2.
CoefName2

#-----
# III. Calculate parameter likelihoods. Do this before W.ParaEst step below where DataSet2 NAs are turned into 0s.

# Create a vector to hold output.
Paralik<- vector()

# Create a matrix to hold output.
ParalikMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(CoefName2), dimnames=list(c(1:length(1)), c(CoefName2))))

# View ParalikMatrix.
ParalikMatrix

# Fill ParalikMatrix manually, b/c loop doesnt work with the subset function!

sub1<- subset(DataSet2, !is.na(Intercept))
nrow(sub1)
ParalikIntercept<- sum(sub1$AIC.c.W)
ParalikMatrix[1]<- signif(ParalikIntercept, digits=12)

sub2<- subset(DataSet2, !is.na(SppCHPE))
nrow(sub2)
ParalikSppCHPE<- sum(sub2$AIC.c.W)
ParalikMatrix[2]<- signif(ParalikSppCHPE, digits=12)

sub3<- subset(DataSet2, !is.na(SppGEPE))
nrow(sub3)
ParalikSppGEPE<- sum(sub3$AIC.c.W)
ParalikMatrix[3]<- signif(ParalikSppGEPE, digits=12)

sub4<- subset(DataSet2, !is.na(Year.22))
nrow(sub4)
ParalikYear.22<- sum(sub4$AIC.c.W)
ParalikMatrix[4]<- signif(ParalikYear.22, digits=12)

sub5<- subset(DataSet2, !is.na(Year.23))
nrow(sub5)
ParalikYear.23<- sum(sub5$AIC.c.W)
ParalikMatrix[5]<- signif(ParalikYear.23, digits=12)

sub6<- subset(DataSet2, !is.na(SppCHPE.Year.22))
nrow(sub6)
ParalikSppCHPE.Year.22<- sum(sub6$AIC.c.W)
ParalikMatrix[6]<- signif(ParalikSppCHPE.Year.22, digits=12)

sub7<- subset(DataSet2, !is.na(SppCHPE.Year.23))
nrow(sub7)
ParalikSppCHPE.Year.23<- sum(sub7$AIC.c.W)
ParalikMatrix[7]<- signif(ParalikSppCHPE.Year.23, digits=12)

sub8<- subset(DataSet2, !is.na(SppGEPE.Year.22))
nrow(sub8)
ParalikSppGEPE.Year.22<- sum(sub8$AIC.c.W)
ParalikMatrix[8]<- signif(ParalikSppGEPE.Year.22, digits=12)

sub9<- subset(DataSet2, !is.na(SppGEPE.Year.23))
nrow(sub9)
ParalikSppGEPE.Year.23<- sum(sub9$AIC.c.W)
ParalikMatrix[9]<- signif(ParalikSppGEPE.Year.23, digits=12)


#-----
# View ParalikMatrix
ParalikMatrix

#-----
# Transpose ParalikMatrix
transpose.ParalikMatrix<- t(ParalikMatrix)

# View transposed ParalikMatrix
transpose.ParalikMatrix

#-----
# Calculate Weighted Parameter Estimates (W.ParaEst).

# Turns NAs into 0s for calculating W.ParaEsts.
DataSet2[is.na(DataSet2)]<- 0

# View DataSet2 to check that NAs have have been replaced.
head(DataSet2)

# Specify a vector for new W.ParaEst.
W.CoefName2<- c("W.Intercept", "W.SppCHPE", "W.SppGEPE", "W.Year.22", "W.Year.23", "W.SppCHPE.Year.22", "W.SppCHPE.Year.23", "W.SppGEPE.Year.22", "W.SppGEPE.Year.23")

# View W.CoefName2.
W.CoefName2

# II. Create matrix to hold outputs for W.ParaEsts.

# Weighted ParaEst Matrix, ncol=8 includes all W.CoefName2 listed above.
W.ParaEstMatrix<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(W.CoefName2), dimnames=list(c(1:nrow(DataSet2)), c(W.CoefName2))))

# View W.ParaMatrix.
head(W.ParaEstMatrix)

# Run W.ParaEst loop to create W.ParaEsts by mulitplying the Para Ests for each model by the AIC.c weight of each model.

for(i in 1:length(CoefName2)) {
  
  W.ParaEst<- DataSet2[CoefName2[i]]*DataSet2["AIC.c.W"]
  W.ParaEstMatrix[,i]<- W.ParaEst
}

# View filled W.ParaEstMatrix.
W.ParaEstMatrix

#-----
# Summed W.ParaEsts.
# Specify a vector for summed W.ParaEsts.
s.W.CoefName2<- c("s.W.Intercept", "s.W.SppCHPE", "s.W.SppGEPE", "s.W.Year.22", "s.W.Year.23", "s.W.SppCHPE.Year.22", "s.W.SppCHPE.Year.23", "s.W.SppGEPE.Year.22", "s.W.SppGEPE.Year.23")

# View summed W.ParaEsts.
s.W.CoefName2

# Create 2 matrices to hold same output for summed W.ParaEsts

# II. 1-Summed W.ParaEsts Matrix, ncol=8 includes all s.W.CoefName2 listed above and nrow(DataSet2) to be used in the W.ParaEst.SE calcs.
s.W.ParaEstMatrix1<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(s.W.CoefName2), dimnames=list(c(1:nrow(DataSet2)), c(s.W.CoefName2))))

# 1-View s.W.ParaEstMatrix
s.W.ParaEstMatrix1

# 1-Run s.W.ParaEst loop. Calculates summed W.ParaEsts by summing all W.ParaEsts.

for(i in 1:length(W.CoefName2)) {
  
  s.W.ParaEst<- sum(W.ParaEstMatrix[W.CoefName2[i]])
  s.W.ParaEst[i]<- s.W.ParaEst
  s.W.ParaEstMatrix1[,i]<- s.W.ParaEst[i]
}

# View filled s.W.ParaEstMatrix1
s.W.ParaEstMatrix1

# III. 2-Summed W.ParaEsts Matrix, ncol=8 includes all s.W.CoefName2 listed above and nrow(1) to be used in final binding of all summed matrices.
s.W.ParaEstMatrix2<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName2), dimnames=list(c(1:length(1)), c(s.W.CoefName2))))

# 2-View s.W.ParaEstMatrix
s.W.ParaEstMatrix2

# 2-Run s.W.ParaEst loop. Calculates summed W.ParaEsts by summing all

for(i in 1:length(W.CoefName2)) {
  
  s.W.ParaEst<- sum(W.ParaEstMatrix[W.CoefName2[i]])
  s.W.ParaEst[i]<- s.W.ParaEst
  s.W.ParaEstMatrix2[,i]<- s.W.ParaEst[i]
}

# 2-View filled s.W.ParaEstMatrix2
s.W.ParaEstMatrix2

# 2-Transpose s.W.ParaEstMatrix2
transpose.s.W.ParaEstMatrix2<- t(s.W.ParaEstMatrix2)

# 2-View transposed s.W.ParaEstMatrix2
transpose.s.W.ParaEstMatrix2

#-----
# Calculate Weighted Parameter Estimate SEs (W.ParaEst.SE).

# Specify a vector for DataSet2 Parameter SEs. These will be used to define the parameter SEs that will used to calculate W.ParaEst.SE that are caculated below in the loop. !!Don't forget to use . instead of : for interactions as this is how they are uploaded in DataSet2.
CoefName3<- c("Intercept.SE", "SppCHPE.SE", "SppGEPE.SE", "Year.22.SE", "Year.23.SE", "SppCHPE.Year.22.SE", "SppCHPE.Year.23.SE", "SppGEPE.Year.22.SE", "SppGEPE.Year.23.SE")

# View CoefName3
CoefName3

# Specify a vector for new W.ParaEst.SE.
W.CoefName3.SE<- c("W.Intercept.SE", "W.SppCHPE.SE", "W.SppGEPE.SE", "W.Year.22.SE", "W.Year.23.SE", "W.SppCHPE.Year.22.SE", "W.SppCHPE.Year.23.SE", "W.SppGEPE.Year.22.SE", "W.SppGEPE.Year.23.SE")

# View W.CoefName2.SE
W.CoefName3.SE

# Create matrix to hold outputs for W.ParaEst.SE.

# II. Weighted ParaEst.SE Matrix, ncol=8 includes all W.CoefName3.SE listed above.
W.ParaEst.SEMatrix<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(W.CoefName3.SE), dimnames=list(c(1:nrow(DataSet2)), c(W.CoefName3.SE))))

# View W.ParaMatrix.
head(W.ParaEst.SEMatrix)

# Run W.ParaEst.SE loop to create W.ParaEstSEs using the following equation W.SE = (AIC.c.W*(sqrt((ParaEst SE^2)+(ParaEst-summedWParaEst)^2)))

for(j in 1:nrow(DataSet2)) {
  for(i in 1:length(CoefName3)) {
    
    if (DataSet2[j,CoefName3[i]] !=0) {W.ParaEst.SEMatrix[j,i]<- DataSet2$AIC.c.W[j]*(sqrt((DataSet2[j,CoefName3[i]])^2+(DataSet2[j,CoefName2[i]]-s.W.ParaEstMatrix1[j,i])^2))} else { 
      W.ParaEst.SEMatrix[j,i]<- 0}
  }
}

# View filled W.ParaEst.SEMatrix.
W.ParaEst.SEMatrix

#-----
# Summed W.ParaEst.SEs (Unconditional SEs)
# Specify a vector for summed W.ParaEsts.SEs.
s.W.CoefName3.SE<- c("s.W.Intercept.SE", "s.W.SppCHPE.SE", "s.W.SppGEPE.SE", "s.W.Year.22.SE", "s.W.Year.23.SE", "s.W.SppCHPE.Year.22.SE", "s.W.SppCHPE.Year.23.SE", "s.W.SppGEPE.Year.22.SE", "s.W.SppGEPE.Year.23.SE")

# View summed W.ParaEst.SEs.
s.W.CoefName3.SE

# Create matrix to hold outputs for summed W.ParaEst.SEs

# III. Summed W.ParaEst.SE Matrix, ncol=8 includes all s.W.CoefNames2.SE listed above
s.W.ParaEst.SEMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName3.SE), dimnames=list(c(1:length(1)), c(s.W.CoefName3.SE))))

# View s.W.ParaEst.SEMatrix
s.W.ParaEst.SEMatrix

#-----
# Run s.W.ParaEst.SE loop. Calculates summed W.ParaEsts by summing all

for(i in 1:length(s.W.CoefName3.SE)) {
  
  s.W.ParaEst.SE<- sum(W.ParaEst.SEMatrix[W.CoefName3.SE[i]])
  s.W.ParaEst.SE[i]<- s.W.ParaEst.SE
  s.W.ParaEst.SEMatrix[,i]<- s.W.ParaEst.SE[i]
}

# View filled s.W.ParaEst.SEMatrix
s.W.ParaEst.SEMatrix

# Transpose s.W.ParaEst.SEMatrix
transpose.s.W.ParaEst.SEMatrix<- t(s.W.ParaEst.SEMatrix)

# View transposed s.W.ParaEst.SEMatrix
transpose.s.W.ParaEst.SEMatrix

#-----
# Unconditional CIs
# Specify a vector for summed W.ParaEst.CIs.
s.W.CoefName3.CI<- c("s.W.Intercept.CI", "s.W.SppCHPE.CI", "s.W.SppGEPE.CI", "s.W.Year.22.CI", "s.W.Year.23.CI", "s.W.SppCHPE.Year.22.CI", "s.W.SppCHPE.Year.23.CI", "s.W.SppGEPE.Year.22.CI", "s.W.SppGEPE.Year.23.CI")

# View summed W.ParaEst.CIs.
s.W.CoefName3.CI

# Create matrix to hold outputs for summed W.ParaEst.CIs

# III. Summed W.ParaEst.CI Matrix, ncol=8 includes all s.W.CoefNames2.CI listed above
s.W.ParaEst.CIMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName3.CI), dimnames=list(c(1:length(1)), c(s.W.CoefName3.CI))))

# View s.W.ParaEst.CIMatrix
s.W.ParaEst.CIMatrix

# Run s.W.ParaEst.CI loop. Calculates summed W.ParaEst.CIs by multiplying 1.96*summed W.ParaEst.Ses

for(i in 1:length(s.W.CoefName3.SE)) {
  
  s.W.ParaEst.CI<- s.W.ParaEst.SEMatrix[s.W.CoefName3.SE[i]]*1.96
  s.W.ParaEst.CIMatrix[,i]<- s.W.ParaEst.CI
}

# View filled s.W.ParaEst.CIMatrix
s.W.ParaEst.CIMatrix

# Transpose s.W.ParaEst.CIMatrix
transpose.s.W.ParaEst.CIMatrix<- t(s.W.ParaEst.CIMatrix)

# View transposed s.W.ParaEst.CIMatrix
transpose.s.W.ParaEst.CIMatrix

#-----
# II. Bind DataSet2, W.ParaEstMatrix, W.ParaEst.SEMatrix

W.AICFinalMatrix<- cbind(DataSet2,W.ParaEstMatrix,W.ParaEst.SEMatrix)

# View W.AICFinaMatrix
W.AICFinalMatrix

# Write W.AICFinaMatrix to a .cvs file
write.table(W.AICFinalMatrix, file="W.AICFinalMatrix.csv", col.names=NA, sep=",")

#-----
# III. Bind ParalikMatrix, s.W.ParaEst.SEMatrix, s.W.ParaEst.CIMatrix
s.W.AICFinalMatrix<- cbind(transpose.ParalikMatrix,transpose.s.W.ParaEstMatrix2,transpose.s.W.ParaEst.SEMatrix,transpose.s.W.ParaEst.CIMatrix)

# View s.W.AICFinalMatrix
s.W.AICFinalMatrix

# Create vector for s.W.AICFinalMatrix column names.
s.W.AICFinalMatrixColName<- c("Paralik","s.W.ParaEst","UncondSE","UncondCI")

s.W.AICFinalMatrixColName

# Write s.W.AICFinaMatrix to a .cvs file. Note will have to move Column names over 1 column to be correct because R is weird.
write.table(s.W.AICFinalMatrix, file="s.W.AICFinalMatrix.csv", col.names=(s.W.AICFinalMatrixColName), sep=",")

#-----
# End code
# You should have now produced the following files with this code:
# 1-AICFinalMatrix.csv; use this file to sort models by deltaAIC.c and in excel calculate manually the cumAIC.c.W values.
# 2-Model Summaries_.doc; this is a word doc that holds the output from all models in the Candidate Set.
# 3-W.AICFinalMatrixWS.csv; this file is used as DataSet2 to run weighted calcs.
# 4-W.AICFinalMatrix.csv; this file contains all Models, Para Ests and associated SEs, AIC.c.Ws for each model, calculated W.Para Ests and associated W.SEs. These are the data that are used to create summed weighted estimates.
# 5-s.W.AICFinalMatrix.csv; this file holds summed weighted estimates including Parameter Likelihoods, s.W.ParaEsts, unconditional SEs, and unconditional CIs.
# !When all these files have been created, cut and paste them into 1 excel file with 4 sheets for each part of the analysis and label this excel file the same name as that for the Model Summaries_.doc file so that things can be referenced easily within the same analysis!

#-----
# Residual plots of most parameterized models.

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see what files are available.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the data to be used.
DataFileName<- "Ads.INS.DataSet.Sub4.csv" #Ads C
DataFileName<- "Ads.INS.DataSet.Sub3.csv" #Ads N

#-----
# Import the data. Be sure to check the structure of the data and that R is reading continuous and categorical variables appropriately.
DataSet<-read.csv(DataFileName)

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

# Change parameter structures if needed.
DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

#-----
mcheck<- function(obj,...) {
  rs<- obj$resid
  fv<- obj$fitted
  par(mfrow=c(1,2))
  plot(fv,rs,xlab="Fitted Values",ylab="Residuals")
  abline(h=0,lty=2)
  qqnorm(rs,xlab="Normal scores",ylab="Ordered residuals",main="")
  par(mfrow=c(1,1))
  invisible(NULL)
}

#-----
# Residual plots of most parameterized model. Carbon.

#-----
DataSet

mod4<- lm(delta.13.C~Spp + Year.2 + Spp:Year.2, data=DataSet)
summary(mod4)

shapiro.test(resid(mod4))

#-----
# jpgs

jpeg(file="ResidPlot delta.13.C~Spp + Year.2 + Spp:Year.2.jpg")

mcheck(mod4)

dev.off()

jpeg(file="HistPlot.Resid delta.13.C~Spp + Year.2 + Spp:Year.2.jpg")

hist(resid(mod4))

dev.off()

# Residual plots of most parameterized model. Nitrogen
mod1<- lm(delta.15.N~Spp + Year.2 + Spp:Year.2, data=DataSet)
summary(mod1)

shapiro.test(resid(mod1))

#-----
# jpgs

jpeg(file="ResidPlot delta.15.N~Spp + Year.2 + Spp:Year.2.jpg")

mcheck(mod1)

dev.off()

jpeg(file="HistPlot.Resid delta.15.N~Spp + Year.2 + Spp:Year.2.jpg")

hist(resid(mod1))

dev.off()

#-----
# Variation in C/N between species analyses. D5 chicks, IDA analysis 2.
# These analyses use a linear mixed effects models because many of the chicks are pairs that come from the same nest,
# so need to control for the fact that their isotope signature is NOT independent.

# LME in R. Use library nlme, load this for each R session. Also the package AICcmodavg produces AIC tables for MEmodels. Mixed effects models are used when the data have some hierarchical form such as in repeated measures, time series and blocked experiments. 

#-----
# Reset R, clear all objects

rm(list=ls()) 

# Packages

install.packages("nlme")
library(nlme)
install.packages("AICcmodavg")
library(AICcmodavg)

#-----
# Load Data

setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

Cks.D5.PAL.INS.DataSet.Sub1<-read.csv("Cks.D5.PAL.INS.DataSet.Sub2.csv") #d5Ch, C & N

head(Cks.D5.PAL.INS.DataSet.Sub1)
nrow(Cks.D5.PAL.INS.DataSet.Sub1)
str(Cks.D5.PAL.INS.DataSet.Sub1)

Cks.D5.PAL.INS.DataSet.Sub1$Year.2<- as.factor(Cks.D5.PAL.INS.DataSet.Sub1$Year.2)

Cks.D5.PAL.INS.DataSet.Sub1$Nest<- as.factor(Cks.D5.PAL.INS.DataSet.Sub1$Nest)

str(Cks.D5.PAL.INS.DataSet.Sub1)

#-----
# Checking lm with no random effect vs lme with random effect for D5 chicks.

mod1<- gls(delta.15.N~Spp + Year.2 + Spp*Year.2, data=Cks.D5.PAL.INS.DataSet.Sub1, method="REML")
mod1
summary(mod1)
AIC.mod1<- AIC(mod1)
AIC.mod1
AIC.c.gls.mod1<- AICc(mod1)
AIC.c.gls.mod1

mod2<- lme(delta.15.N~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="REML")
summary(mod2)
AIC.mod2<- AIC(mod2)
AIC.mod2
AIC.c.mod2<- AICc(mod2)
AIC.c.mod2

mod3<- gls(delta.13.C~Spp + Year.2 + Spp*Year.2, data=Cks.D5.PAL.INS.DataSet.Sub1, method="REML")
mod3
summary(mod3)
AIC.mod3<- AIC(mod3)
AIC.mod3
AIC.c.lm.mod3<- AICc(mod3)
AIC.c.lm.mod3

mod4<- lme(delta.13.C~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="REML")
summary(mod4)
AIC.mod4<- AIC(mod4)
AIC.mod4
AIC.c.mod4<- AICc(mod4)
AIC.c.mod4

#-----
# Carbon
str(Cks.D5.PAL.INS.DataSet.Sub1)

Cand.models<- list()

Cand.models[[1]]<- lme(delta.13.C~1, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[1]])
resid.mod1<- resid(Cand.models[[1]])
resid.mod1.2<- resid(Cand.models[[1]])^2
resid.mod1.2.sum<- sum(resid.mod1.2)
mod1.r2<- 1-(resid.mod1.2.sum/resid.mod1.2.sum)

Cand.models[[2]]<- lme(delta.13.C~Spp, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[2]])
resid.mod2<- resid(Cand.models[[2]])
resid.mod2.2<- resid(Cand.models[[2]])^2
resid.mod2.2.sum<- sum(resid.mod2.2)
mod2.r2<- 1-(resid.mod2.2.sum/resid.mod1.2.sum)

Cand.models[[3]]<- lme(delta.13.C~Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[3]])
resid.mod3<- resid(Cand.models[[3]])
resid.mod3.2<- resid(Cand.models[[3]])^2
resid.mod3.2.sum<- sum(resid.mod3.2)
mod3.r2<- 1-(resid.mod3.2.sum/resid.mod1.2.sum)

Cand.models[[4]]<- lme(delta.13.C~Spp + Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[4]])
resid.mod4<- resid(Cand.models[[4]])
resid.mod4.2<- resid(Cand.models[[4]])^2
resid.mod4.2.sum<- sum(resid.mod4.2)
mod4.r2<- 1-(resid.mod4.2.sum/resid.mod1.2.sum)

Cand.models[[5]]<- lme(delta.13.C~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[5]])
resid.mod5<- resid(Cand.models[[5]])
resid.mod5.2<- resid(Cand.models[[5]])^2
resid.mod5.2.sum<- sum(resid.mod5.2)
mod5.r2<- 1-(resid.mod5.2.sum/resid.mod1.2.sum)

Modnames<- paste("mod", 1:length(Cand.models), sep="")

AICtable<- aictab(cand.set = Cand.models, modnames = Modnames, sort = FALSE)
AICtable

AICtable.2<- print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),digits = 5, LL = TRUE)

Cks.D5.PAL.INS.Resid.R2<- cbind(resid.mod1, resid.mod1.2, resid.mod1.2.sum, mod1.r2, resid.mod2, resid.mod2.2, resid.mod2.2.sum, mod2.r2, resid.mod3, resid.mod3.2, resid.mod3.2.sum, mod3.r2, resid.mod4, resid.mod4.2, resid.mod4.2.sum, mod4.r2, resid.mod5, resid.mod5.2, resid.mod5.2.sum, mod5.r2)

Cks.D5.PAL.INS.Resid.R2

write.table(Cks.D5.PAL.INS.Resid.R2, file="Cks.D5.PAL.INS.Resid.C.R2.csv", col.names=NA, sep=",")

#-----
# Nitrogen
str(Cks.D5.PAL.INS.DataSet.Sub1)

Cand.models<- list()

Cand.models[[1]]<- lme(delta.15.N~1, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[1]])
resid.mod1<- resid(Cand.models[[1]])
resid.mod1.2<- resid(Cand.models[[1]])^2
resid.mod1.2.sum<- sum(resid.mod1.2)
mod1.r2<- 1-(resid.mod1.2.sum/resid.mod1.2.sum)

Cand.models[[2]]<- lme(delta.15.N~Spp, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[2]])
resid.mod2<- resid(Cand.models[[2]])
resid.mod2.2<- resid(Cand.models[[2]])^2
resid.mod2.2.sum<- sum(resid.mod2.2)
mod2.r2<- 1-(resid.mod2.2.sum/resid.mod1.2.sum)

Cand.models[[3]]<- lme(delta.15.N~Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[3]])
resid.mod3<- resid(Cand.models[[3]])
resid.mod3.2<- resid(Cand.models[[3]])^2
resid.mod3.2.sum<- sum(resid.mod3.2)
mod3.r2<- 1-(resid.mod3.2.sum/resid.mod1.2.sum)

Cand.models[[4]]<- lme(delta.15.N~Spp + Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[4]])
resid.mod4<- resid(Cand.models[[4]])
resid.mod4.2<- resid(Cand.models[[4]])^2
resid.mod4.2.sum<- sum(resid.mod4.2)
mod4.r2<- 1-(resid.mod4.2.sum/resid.mod1.2.sum)

Cand.models[[5]]<- lme(delta.15.N~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D5.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[5]])
resid.mod5<- resid(Cand.models[[5]])
resid.mod5.2<- resid(Cand.models[[5]])^2
resid.mod5.2.sum<- sum(resid.mod5.2)
mod5.r2<- 1-(resid.mod5.2.sum/resid.mod1.2.sum)

Modnames<- paste("mod", 1:length(Cand.models), sep="")

AICtable<- aictab(cand.set = Cand.models, modnames = Modnames, sort = FALSE)
AICtable

AICtable.2<- print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),digits = 5, LL = TRUE)

Cks.D5.PAL.INS.Resid.R2<- cbind(resid.mod1, resid.mod1.2, resid.mod1.2.sum, mod1.r2, resid.mod2, resid.mod2.2, resid.mod2.2.sum, mod2.r2, resid.mod3, resid.mod3.2, resid.mod3.2.sum, mod3.r2, resid.mod4, resid.mod4.2, resid.mod4.2.sum, mod4.r2, resid.mod5, resid.mod5.2, resid.mod5.2.sum, mod5.r2)

Cks.D5.PAL.INS.Resid.R2

write.table(Cks.D5.PAL.INS.Resid.R2, file="Cks.D5.PAL.INS.Resid.N.R2.csv", col.names=NA, sep=",")

#-----
# Model averaging

modellist<-list(Cand.models[[1]],Cand.models[[2]],Cand.models[[3]],Cand.models[[4]],Cand.models[[5]])

modnames<- c("1","2","3","4","5")

par.names<-c("Intercept", "SppCHPE", "SppGEPE", "Year.22", "Year.23", "SppCHPE:Year.22", "SppGEPE:Year.22", "SppCHPE:Year.23", "SppGEPE:Year.23")

par.names.se<-c("Intercept.SE", "SppCHPE.SE", "SppGEPE.SE", "Year.22.SE", "Year.23.SE", "SppCHPE:Year.22.SE", "SppGEPE:Year.22.SE", "SppCHPE:Year.23.SE", "SppGEPE:Year.23.SE")

CoefMatrix<- matrix(NA, nrow=length(modnames), ncol=9, dimnames = list(modnames, par.names))
CoefMatrix

SEMatrix<-matrix(NA, nrow=length(modnames), ncol=9, dimnames = list(modnames, par.names.se))
SEMatrix

CoefMatrix[1,1] <- summary(modellist[[1]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[2,1] <- summary(modellist[[2]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[2,2] <- summary(modellist[[2]])$tTable[c("SppCHPE"),c("Value")]
CoefMatrix[2,3] <- summary(modellist[[2]])$tTable[c("SppGEPE"),c("Value")]
CoefMatrix[3,1] <- summary(modellist[[3]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[3,4] <- summary(modellist[[3]])$tTable[c("Year.22"),c("Value")]
CoefMatrix[3,5] <- summary(modellist[[3]])$tTable[c("Year.23"),c("Value")]
CoefMatrix[4,1] <- summary(modellist[[4]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[4,2] <- summary(modellist[[4]])$tTable[c("SppCHPE"),c("Value")]
CoefMatrix[4,3] <- summary(modellist[[4]])$tTable[c("SppGEPE"),c("Value")]
CoefMatrix[4,4] <- summary(modellist[[4]])$tTable[c("Year.22"),c("Value")]
CoefMatrix[4,5] <- summary(modellist[[4]])$tTable[c("Year.23"),c("Value")]
CoefMatrix[5,1] <- summary(modellist[[5]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[5,2] <- summary(modellist[[5]])$tTable[c("SppCHPE"),c("Value")]
CoefMatrix[5,3] <- summary(modellist[[5]])$tTable[c("SppGEPE"),c("Value")]
CoefMatrix[5,4] <- summary(modellist[[5]])$tTable[c("Year.22"),c("Value")]
CoefMatrix[5,5] <- summary(modellist[[5]])$tTable[c("Year.23"),c("Value")]
CoefMatrix[5,6] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.22"),c("Value")]
CoefMatrix[5,7] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.22"),c("Value")]
CoefMatrix[5,8] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.23"),c("Value")]
CoefMatrix[5,9] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.23"),c("Value")]

SEMatrix[1,1] <- summary(modellist[[1]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[2,1] <- summary(modellist[[2]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[2,2] <- summary(modellist[[2]])$tTable[c("SppCHPE"),c("Std.Error")]
SEMatrix[2,3] <- summary(modellist[[2]])$tTable[c("SppGEPE"),c("Std.Error")]
SEMatrix[3,1] <- summary(modellist[[3]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[3,4] <- summary(modellist[[3]])$tTable[c("Year.22"),c("Std.Error")]
SEMatrix[3,5] <- summary(modellist[[3]])$tTable[c("Year.23"),c("Std.Error")]
SEMatrix[4,1] <- summary(modellist[[4]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[4,2] <- summary(modellist[[4]])$tTable[c("SppCHPE"),c("Std.Error")]
SEMatrix[4,3] <- summary(modellist[[4]])$tTable[c("SppGEPE"),c("Std.Error")]
SEMatrix[4,4] <- summary(modellist[[4]])$tTable[c("Year.22"),c("Std.Error")]
SEMatrix[4,5] <- summary(modellist[[4]])$tTable[c("Year.23"),c("Std.Error")]
SEMatrix[5,1] <- summary(modellist[[5]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[5,2] <- summary(modellist[[5]])$tTable[c("SppCHPE"),c("Std.Error")]
SEMatrix[5,3] <- summary(modellist[[5]])$tTable[c("SppGEPE"),c("Std.Error")]
SEMatrix[5,4] <- summary(modellist[[5]])$tTable[c("Year.22"),c("Std.Error")]
SEMatrix[5,5] <- summary(modellist[[5]])$tTable[c("Year.23"),c("Std.Error")]
SEMatrix[5,6] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.22"),c("Std.Error")]
SEMatrix[5,7] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.22"),c("Std.Error")]
SEMatrix[5,8] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.23"),c("Std.Error")]
SEMatrix[5,9] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.23"),c("Std.Error")]

CoefMatrix.filled<- CoefMatrix
CoefMatrix.filled

SEMatrix.filled<- SEMatrix
SEMatrix.filled

#-----
# Parameter likelihoods

Paralik<- vector()

ParalikMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names)))
ParalikMatrix

AICtable

ParalikMatrix[1]<-sum(AICtable[1:5,c("AICcWt")]) # Intercept
ParalikMatrix[2]<-sum(AICtable[c(2,4,5),c("AICcWt")]) # SppCHPE
ParalikMatrix[3]<-sum(AICtable[c(2,4,5),c("AICcWt")]) # SppGEPE
ParalikMatrix[4]<-sum(AICtable[c(3,4,5),c("AICcWt")]) # Year.22
ParalikMatrix[5]<-sum(AICtable[c(3,4,5),c("AICcWt")]) # Year.23
ParalikMatrix[6]<-sum(AICtable[c(5),c("AICcWt")]) # SppCHPE:Year.22
ParalikMatrix[7]<-sum(AICtable[c(5),c("AICcWt")]) # SppGEPE:Year.22
ParalikMatrix[8]<-sum(AICtable[c(5),c("AICcWt")]) # SppCHPE:Year.23
ParalikMatrix[9]<-sum(AICtable[c(5),c("AICcWt")]) # SppGEPE:Year.23

ParalikMatrix.filled<- ParalikMatrix
ParalikMatrix.filled

t(ParalikMatrix.filled)

#-----
# Weighted parameter estimates

AICtable

CoefMatrix.filled[is.na(CoefMatrix.filled)] <- 0
CoefMatrix.filled

s.W.ParaEstMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names)))
s.W.ParaEstMatrix

s.W.ParaEstMatrix[1]<-sum(CoefMatrix.filled[,c("Intercept")]*AICtable$AICcWt)
s.W.ParaEstMatrix[2]<-sum(CoefMatrix.filled[,c("SppCHPE")]*AICtable$AICcWt)
s.W.ParaEstMatrix[3]<-sum(CoefMatrix.filled[,c("SppGEPE")]*AICtable$AICcWt)
s.W.ParaEstMatrix[4]<-sum(CoefMatrix.filled[,c("Year.22")]*AICtable$AICcWt)
s.W.ParaEstMatrix[5]<-sum(CoefMatrix.filled[,c("Year.23")]*AICtable$AICcWt)
s.W.ParaEstMatrix[6]<-sum(CoefMatrix.filled[,c("SppCHPE:Year.22")]*AICtable$AICcWt)
s.W.ParaEstMatrix[7]<-sum(CoefMatrix.filled[,c("SppGEPE:Year.22")]*AICtable$AICcWt)
s.W.ParaEstMatrix[8]<-sum(CoefMatrix.filled[,c("SppCHPE:Year.23")]*AICtable$AICcWt)
s.W.ParaEstMatrix[9]<-sum(CoefMatrix.filled[,c("SppGEPE:Year.23")]*AICtable$AICcWt)

s.W.ParaEstMatrix.filled<- s.W.ParaEstMatrix
s.W.ParaEstMatrix.filled

t(s.W.ParaEstMatrix.filled)

#-----
# Unconditional SEs

AICtable

uncondSEMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names.se)))
uncondSEMatrix

CoefMatrix.filled<- CoefMatrix
CoefMatrix.filled

CoefMatrix.filled[is.na(CoefMatrix.filled)] <- 0
CoefMatrix.filled

SEMatrix.filled<- SEMatrix
SEMatrix.filled

SEMatrix.filled[is.na(SEMatrix.filled)] <- 0
SEMatrix.filled

uncondSEMatrix[1]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"Intercept.SE"]^2)+(CoefMatrix.filled[,c("Intercept")]-s.W.ParaEstMatrix.filled[1])^2)))
uncondSEMatrix[2]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppCHPE.SE"]^2)+(CoefMatrix.filled[,c("SppCHPE")]-s.W.ParaEstMatrix.filled[2])^2)))
uncondSEMatrix[3]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppGEPE.SE"]^2)+(CoefMatrix.filled[,c("SppGEPE")]-s.W.ParaEstMatrix.filled[3])^2)))
uncondSEMatrix[4]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"Year.22.SE"]^2)+(CoefMatrix.filled[,c("Year.22")]-s.W.ParaEstMatrix.filled[4])^2)))
uncondSEMatrix[5]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"Year.23.SE"]^2)+(CoefMatrix.filled[,c("Year.23")]-s.W.ParaEstMatrix.filled[5])^2)))
uncondSEMatrix[6]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppCHPE:Year.22.SE"]^2)+(CoefMatrix.filled[,c("SppCHPE:Year.22")]-s.W.ParaEstMatrix.filled[6])^2)))
uncondSEMatrix[7]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppGEPE:Year.22.SE"]^2)+(CoefMatrix.filled[,c("SppGEPE:Year.22")]-s.W.ParaEstMatrix.filled[7])^2)))
uncondSEMatrix[8]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppCHPE:Year.23.SE"]^2)+(CoefMatrix.filled[,c("SppCHPE:Year.23")]-s.W.ParaEstMatrix.filled[8])^2)))
uncondSEMatrix[9]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppGEPE:Year.23.SE"]^2)+(CoefMatrix.filled[,c("SppGEPE:Year.23")]-s.W.ParaEstMatrix.filled[9])^2)))

uncondSEMatrix.filled<- uncondSEMatrix
uncondSEMatrix.filled

t(uncondSEMatrix.filled)

#-----
# Unconditional CIs

par.names.ci<-c("Intercept.CI", "SppCHPE.CI", "SppGEPE.CI", "Year.22.CI", "Year.23.CI", "SppCHPE:Year.22.CI", "SppGEPE:Year.22.CI", "SppCHPE:Year.23.CI", "SppGEPE:Year.23.CI")

uncondCIMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names.ci)))
uncondCIMatrix

uncondCIMatrix[1]<- uncondSEMatrix.filled[1]*1.96
uncondCIMatrix[2]<- uncondSEMatrix.filled[2]*1.96
uncondCIMatrix[3]<- uncondSEMatrix.filled[3]*1.96
uncondCIMatrix[4]<- uncondSEMatrix.filled[4]*1.96
uncondCIMatrix[5]<- uncondSEMatrix.filled[5]*1.96
uncondCIMatrix[6]<- uncondSEMatrix.filled[6]*1.96
uncondCIMatrix[7]<- uncondSEMatrix.filled[7]*1.96
uncondCIMatrix[8]<- uncondSEMatrix.filled[8]*1.96
uncondCIMatrix[9]<- uncondSEMatrix.filled[9]*1.96

uncondCIMatrix.filled<- uncondCIMatrix
uncondCIMatrix

t(uncondCIMatrix)

#-----
# Summary

AICtable

AICtable.2<- print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),digits = 5, LL = TRUE)

CoefMatrix

SEMatrix

t(ParalikMatrix.filled)

t(s.W.ParaEstMatrix.filled)

t(uncondSEMatrix.filled)

t(uncondCIMatrix)

#-----
# Variation in C/N between species analyses. D15 chicks, IDA analysis 3.
# These analyses use a linear mixed effects models because many of the chicks are pairs that come from the same nest,
# so need to control for the fact that their isotope signature is NOT independent.

# LME in R. Use library nlme, load this for each R session. Also the package AICcmodavg produces AIC tables for MEmodels. Mixed effects models are used when the data have some hierarchical form such as in repeated measures, time series and blocked experiments. 

#-----
# Reset R, clear all objects

rm(list=ls()) 

# Packages

install.packages("nlme")
library(nlme)
install.packages("AICcmodavg")
library(AICcmodavg)

# Load Data

setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

Cks.D15.PAL.INS.DataSet.Sub1<-read.csv("Cks.D15.PAL.INS.DataSet.Sub2.csv")

head(Cks.D15.PAL.INS.DataSet.Sub1)
nrow(Cks.D15.PAL.INS.DataSet.Sub1)
str(Cks.D15.PAL.INS.DataSet.Sub1)

Cks.D15.PAL.INS.DataSet.Sub1$Year.2<- as.factor(Cks.D15.PAL.INS.DataSet.Sub1$Year.2)

Cks.D15.PAL.INS.DataSet.Sub1$Nest<- as.factor(Cks.D15.PAL.INS.DataSet.Sub1$Nest)

str(Cks.D15.PAL.INS.DataSet.Sub1)

#-----
# Checking lm with no random effect vs lme with random effect for D15 chicks.

mod1<- gls(delta.15.N~Spp + Year.2 + Spp*Year.2, data=Cks.D15.PAL.INS.DataSet.Sub1, method="REML")
mod1
summary(mod1)
AIC.mod1<- AIC(mod1)
AIC.mod1
AIC.c.lm.mod1<- AICc(mod1)
AIC.c.lm.mod1

mod2<- lme(delta.15.N~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="REML")
summary(mod2)
AIC.mod2<- AIC(mod2)
AIC.mod2
AIC.c.mod2<- AICc(mod2)
AIC.c.mod2

mod3<- gls(delta.13.C~Spp + Year.2 + Spp*Year.2, data=Cks.D15.PAL.INS.DataSet.Sub1, method="REML")
mod3
summary(mod3)
AIC.mod3<- AIC(mod3)
AIC.mod3
AIC.c.lm.mod3<- AICc(mod3)
AIC.c.lm.mod3

mod4<- lme(delta.13.C~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="REML")
summary(mod4)
AIC.mod4<- AIC(mod4)
AIC.mod4
AIC.c.mod4<- AICc(mod4)
AIC.c.mod4

#-----
# Carbon
str(Cks.D15.PAL.INS.DataSet.Sub1)

Cand.models<- list()

Cand.models[[1]]<- lme(delta.13.C~1, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[1]])
resid.mod1<- resid(Cand.models[[1]])
resid.mod1.2<- resid(Cand.models[[1]])^2
resid.mod1.2.sum<- sum(resid.mod1.2)
mod1.r2<- 1-(resid.mod1.2.sum/resid.mod1.2.sum)

Cand.models[[2]]<- lme(delta.13.C~Spp, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[2]])
resid.mod2<- resid(Cand.models[[2]])
resid.mod2.2<- resid(Cand.models[[2]])^2
resid.mod2.2.sum<- sum(resid.mod2.2)
mod2.r2<- 1-(resid.mod2.2.sum/resid.mod1.2.sum)

Cand.models[[3]]<- lme(delta.13.C~Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[3]])
resid.mod3<- resid(Cand.models[[3]])
resid.mod3.2<- resid(Cand.models[[3]])^2
resid.mod3.2.sum<- sum(resid.mod3.2)
mod3.r2<- 1-(resid.mod3.2.sum/resid.mod1.2.sum)

Cand.models[[4]]<- lme(delta.13.C~Spp + Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[4]])
resid.mod4<- resid(Cand.models[[4]])
resid.mod4.2<- resid(Cand.models[[4]])^2
resid.mod4.2.sum<- sum(resid.mod4.2)
mod4.r2<- 1-(resid.mod4.2.sum/resid.mod1.2.sum)

Cand.models[[5]]<- lme(delta.13.C~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[5]])
resid.mod5<- resid(Cand.models[[5]])
resid.mod5.2<- resid(Cand.models[[5]])^2
resid.mod5.2.sum<- sum(resid.mod5.2)
mod5.r2<- 1-(resid.mod5.2.sum/resid.mod1.2.sum)

Modnames<- paste("mod", 1:length(Cand.models), sep="")

AICtable<- aictab(cand.set = Cand.models, modnames = Modnames, sort = FALSE)
AICtable

AICtable.2<- print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),digits = 5, LL = TRUE)

Cks.D5.PAL.INS.Resid.R2<- cbind(resid.mod1, resid.mod1.2, resid.mod1.2.sum, mod1.r2, resid.mod2, resid.mod2.2, resid.mod2.2.sum, mod2.r2, resid.mod3, resid.mod3.2, resid.mod3.2.sum, mod3.r2, resid.mod4, resid.mod4.2, resid.mod4.2.sum, mod4.r2, resid.mod5, resid.mod5.2, resid.mod5.2.sum, mod5.r2)

Cks.D5.PAL.INS.Resid.R2

write.table(Cks.D5.PAL.INS.Resid.R2, file="Cks.D15.PAL.INS.Resid.R2.C.csv", col.names=NA, sep=",")

#-----
# Nitrogen
str(Cks.D15.PAL.INS.DataSet.Sub1)

Cand.models<- list()

Cand.models[[1]]<- lme(delta.15.N~1, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[1]])
resid.mod1<- resid(Cand.models[[1]])
resid.mod1.2<- resid(Cand.models[[1]])^2
resid.mod1.2.sum<- sum(resid.mod1.2)
mod1.r2<- 1-(resid.mod1.2.sum/resid.mod1.2.sum)

Cand.models[[2]]<- lme(delta.15.N~Spp, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[2]])
resid.mod2<- resid(Cand.models[[2]])
resid.mod2.2<- resid(Cand.models[[2]])^2
resid.mod2.2.sum<- sum(resid.mod2.2)
mod2.r2<- 1-(resid.mod2.2.sum/resid.mod1.2.sum)

Cand.models[[3]]<- lme(delta.15.N~Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[3]])
resid.mod3<- resid(Cand.models[[3]])
resid.mod3.2<- resid(Cand.models[[3]])^2
resid.mod3.2.sum<- sum(resid.mod3.2)
mod3.r2<- 1-(resid.mod3.2.sum/resid.mod1.2.sum)

Cand.models[[4]]<- lme(delta.15.N~Spp + Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[4]])
resid.mod4<- resid(Cand.models[[4]])
resid.mod4.2<- resid(Cand.models[[4]])^2
resid.mod4.2.sum<- sum(resid.mod4.2)
mod4.r2<- 1-(resid.mod4.2.sum/resid.mod1.2.sum)

Cand.models[[5]]<- lme(delta.15.N~Spp + Year.2 + Spp*Year.2, random=~1|Nest, data=Cks.D15.PAL.INS.DataSet.Sub1, method="ML")
summary(Cand.models[[5]])
resid.mod5<- resid(Cand.models[[5]])
resid.mod5.2<- resid(Cand.models[[5]])^2
resid.mod5.2.sum<- sum(resid.mod5.2)
mod5.r2<- 1-(resid.mod5.2.sum/resid.mod1.2.sum)

Modnames<- paste("mod", 1:length(Cand.models), sep="")

AICtable<- aictab(cand.set = Cand.models, modnames = Modnames, sort = FALSE)
AICtable

AICtable.2<- print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),digits = 5, LL = TRUE)

Cks.D5.PAL.INS.Resid.R2<- cbind(resid.mod1, resid.mod1.2, resid.mod1.2.sum, mod1.r2, resid.mod2, resid.mod2.2, resid.mod2.2.sum, mod2.r2, resid.mod3, resid.mod3.2, resid.mod3.2.sum, mod3.r2, resid.mod4, resid.mod4.2, resid.mod4.2.sum, mod4.r2, resid.mod5, resid.mod5.2, resid.mod5.2.sum, mod5.r2)

Cks.D5.PAL.INS.Resid.R2

write.table(Cks.D5.PAL.INS.Resid.R2, file="Cks.D15.PAL.INS.Resid.R2.N.csv", col.names=NA, sep=",")

#-----
# Model averaging

modellist<-list(Cand.models[[1]],Cand.models[[2]],Cand.models[[3]],Cand.models[[4]],Cand.models[[5]])

modnames<- c("1","2","3","4","5")

par.names<-c("Intercept", "SppCHPE", "SppGEPE", "Year.22", "Year.23", "SppCHPE:Year.22", "SppGEPE:Year.22", "SppCHPE:Year.23", "SppGEPE:Year.23")

par.names.se<-c("Intercept.SE", "SppCHPE.SE", "SppGEPE.SE", "Year.22.SE", "Year.23.SE", "SppCHPE:Year.22.SE", "SppGEPE:Year.22.SE", "SppCHPE:Year.23.SE", "SppGEPE:Year.23.SE")

CoefMatrix<- matrix(NA, nrow=length(modnames), ncol=9, dimnames = list(modnames, par.names))
CoefMatrix

SEMatrix<-matrix(NA, nrow=length(modnames), ncol=9, dimnames = list(modnames, par.names.se))
SEMatrix

CoefMatrix[1,1] <- summary(modellist[[1]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[2,1] <- summary(modellist[[2]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[2,2] <- summary(modellist[[2]])$tTable[c("SppCHPE"),c("Value")]
CoefMatrix[2,3] <- summary(modellist[[2]])$tTable[c("SppGEPE"),c("Value")]
CoefMatrix[3,1] <- summary(modellist[[3]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[3,4] <- summary(modellist[[3]])$tTable[c("Year.22"),c("Value")]
CoefMatrix[3,5] <- summary(modellist[[3]])$tTable[c("Year.23"),c("Value")]
CoefMatrix[4,1] <- summary(modellist[[4]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[4,2] <- summary(modellist[[4]])$tTable[c("SppCHPE"),c("Value")]
CoefMatrix[4,3] <- summary(modellist[[4]])$tTable[c("SppGEPE"),c("Value")]
CoefMatrix[4,4] <- summary(modellist[[4]])$tTable[c("Year.22"),c("Value")]
CoefMatrix[4,5] <- summary(modellist[[4]])$tTable[c("Year.23"),c("Value")]
CoefMatrix[5,1] <- summary(modellist[[5]])$tTable[c("(Intercept)"),c("Value")]
CoefMatrix[5,2] <- summary(modellist[[5]])$tTable[c("SppCHPE"),c("Value")]
CoefMatrix[5,3] <- summary(modellist[[5]])$tTable[c("SppGEPE"),c("Value")]
CoefMatrix[5,4] <- summary(modellist[[5]])$tTable[c("Year.22"),c("Value")]
CoefMatrix[5,5] <- summary(modellist[[5]])$tTable[c("Year.23"),c("Value")]
CoefMatrix[5,6] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.22"),c("Value")]
CoefMatrix[5,7] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.22"),c("Value")]
CoefMatrix[5,8] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.23"),c("Value")]
CoefMatrix[5,9] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.23"),c("Value")]

SEMatrix[1,1] <- summary(modellist[[1]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[2,1] <- summary(modellist[[2]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[2,2] <- summary(modellist[[2]])$tTable[c("SppCHPE"),c("Std.Error")]
SEMatrix[2,3] <- summary(modellist[[2]])$tTable[c("SppGEPE"),c("Std.Error")]
SEMatrix[3,1] <- summary(modellist[[3]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[3,4] <- summary(modellist[[3]])$tTable[c("Year.22"),c("Std.Error")]
SEMatrix[3,5] <- summary(modellist[[3]])$tTable[c("Year.23"),c("Std.Error")]
SEMatrix[4,1] <- summary(modellist[[4]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[4,2] <- summary(modellist[[4]])$tTable[c("SppCHPE"),c("Std.Error")]
SEMatrix[4,3] <- summary(modellist[[4]])$tTable[c("SppGEPE"),c("Std.Error")]
SEMatrix[4,4] <- summary(modellist[[4]])$tTable[c("Year.22"),c("Std.Error")]
SEMatrix[4,5] <- summary(modellist[[4]])$tTable[c("Year.23"),c("Std.Error")]
SEMatrix[5,1] <- summary(modellist[[5]])$tTable[c("(Intercept)"),c("Std.Error")]
SEMatrix[5,2] <- summary(modellist[[5]])$tTable[c("SppCHPE"),c("Std.Error")]
SEMatrix[5,3] <- summary(modellist[[5]])$tTable[c("SppGEPE"),c("Std.Error")]
SEMatrix[5,4] <- summary(modellist[[5]])$tTable[c("Year.22"),c("Std.Error")]
SEMatrix[5,5] <- summary(modellist[[5]])$tTable[c("Year.23"),c("Std.Error")]
SEMatrix[5,6] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.22"),c("Std.Error")]
SEMatrix[5,7] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.22"),c("Std.Error")]
SEMatrix[5,8] <- summary(modellist[[5]])$tTable[c("SppCHPE:Year.23"),c("Std.Error")]
SEMatrix[5,9] <- summary(modellist[[5]])$tTable[c("SppGEPE:Year.23"),c("Std.Error")]

CoefMatrix.filled<- CoefMatrix
CoefMatrix.filled

SEMatrix.filled<- SEMatrix
SEMatrix.filled

#-----
# Parameter likelihoods

Paralik<- vector()

ParalikMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names)))
ParalikMatrix

AICtable

ParalikMatrix[1]<-sum(AICtable[1:5,c("AICcWt")]) # Intercept
ParalikMatrix[2]<-sum(AICtable[c(2,4,5),c("AICcWt")]) # SppCHPE
ParalikMatrix[3]<-sum(AICtable[c(2,4,5),c("AICcWt")]) # SppGEPE
ParalikMatrix[4]<-sum(AICtable[c(3,4,5),c("AICcWt")]) # Year.22
ParalikMatrix[5]<-sum(AICtable[c(3,4,5),c("AICcWt")]) # Year.23
ParalikMatrix[6]<-sum(AICtable[c(5),c("AICcWt")]) # SppCHPE:Year.22
ParalikMatrix[7]<-sum(AICtable[c(5),c("AICcWt")]) # SppGEPE:Year.22
ParalikMatrix[8]<-sum(AICtable[c(5),c("AICcWt")]) # SppCHPE:Year.23
ParalikMatrix[9]<-sum(AICtable[c(5),c("AICcWt")]) # SppGEPE:Year.23

ParalikMatrix.filled<- ParalikMatrix
ParalikMatrix.filled

t(ParalikMatrix.filled)

#-----
# Weighted parameter estimates

AICtable

CoefMatrix.filled[is.na(CoefMatrix.filled)] <- 0
CoefMatrix.filled

s.W.ParaEstMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names)))
s.W.ParaEstMatrix

s.W.ParaEstMatrix[1]<-sum(CoefMatrix.filled[,c("Intercept")]*AICtable$AICcWt)
s.W.ParaEstMatrix[2]<-sum(CoefMatrix.filled[,c("SppCHPE")]*AICtable$AICcWt)
s.W.ParaEstMatrix[3]<-sum(CoefMatrix.filled[,c("SppGEPE")]*AICtable$AICcWt)
s.W.ParaEstMatrix[4]<-sum(CoefMatrix.filled[,c("Year.22")]*AICtable$AICcWt)
s.W.ParaEstMatrix[5]<-sum(CoefMatrix.filled[,c("Year.23")]*AICtable$AICcWt)
s.W.ParaEstMatrix[6]<-sum(CoefMatrix.filled[,c("SppCHPE:Year.22")]*AICtable$AICcWt)
s.W.ParaEstMatrix[7]<-sum(CoefMatrix.filled[,c("SppGEPE:Year.22")]*AICtable$AICcWt)
s.W.ParaEstMatrix[8]<-sum(CoefMatrix.filled[,c("SppCHPE:Year.23")]*AICtable$AICcWt)
s.W.ParaEstMatrix[9]<-sum(CoefMatrix.filled[,c("SppGEPE:Year.23")]*AICtable$AICcWt)

s.W.ParaEstMatrix.filled<- s.W.ParaEstMatrix
s.W.ParaEstMatrix.filled

t(s.W.ParaEstMatrix.filled)

#-----
# Unconditional SEs

AICtable

uncondSEMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names.se)))
uncondSEMatrix

CoefMatrix.filled<- CoefMatrix
CoefMatrix.filled

CoefMatrix.filled[is.na(CoefMatrix.filled)] <- 0
CoefMatrix.filled

SEMatrix.filled<- SEMatrix
SEMatrix.filled

SEMatrix.filled[is.na(SEMatrix.filled)] <- 0
SEMatrix.filled

uncondSEMatrix[1]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"Intercept.SE"]^2)+(CoefMatrix.filled[,c("Intercept")]-s.W.ParaEstMatrix.filled[1])^2)))
uncondSEMatrix[2]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppCHPE.SE"]^2)+(CoefMatrix.filled[,c("SppCHPE")]-s.W.ParaEstMatrix.filled[2])^2)))
uncondSEMatrix[3]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppGEPE.SE"]^2)+(CoefMatrix.filled[,c("SppGEPE")]-s.W.ParaEstMatrix.filled[3])^2)))
uncondSEMatrix[4]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"Year.22.SE"]^2)+(CoefMatrix.filled[,c("Year.22")]-s.W.ParaEstMatrix.filled[4])^2)))
uncondSEMatrix[5]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"Year.23.SE"]^2)+(CoefMatrix.filled[,c("Year.23")]-s.W.ParaEstMatrix.filled[5])^2)))
uncondSEMatrix[6]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppCHPE:Year.22.SE"]^2)+(CoefMatrix.filled[,c("SppCHPE:Year.22")]-s.W.ParaEstMatrix.filled[6])^2)))
uncondSEMatrix[7]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppGEPE:Year.22.SE"]^2)+(CoefMatrix.filled[,c("SppGEPE:Year.22")]-s.W.ParaEstMatrix.filled[7])^2)))
uncondSEMatrix[8]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppCHPE:Year.23.SE"]^2)+(CoefMatrix.filled[,c("SppCHPE:Year.23")]-s.W.ParaEstMatrix.filled[8])^2)))
uncondSEMatrix[9]<-sum(AICtable$AICcWt*(sqrt((SEMatrix.filled[,"SppGEPE:Year.23.SE"]^2)+(CoefMatrix.filled[,c("SppGEPE:Year.23")]-s.W.ParaEstMatrix.filled[9])^2)))

uncondSEMatrix.filled<- uncondSEMatrix
uncondSEMatrix.filled

t(uncondSEMatrix.filled)

#-----
# Unconditional CIs

par.names.ci<-c("Intercept.CI", "SppCHPE.CI", "SppGEPE.CI", "Year.22.CI", "Year.23.CI", "SppCHPE:Year.22.CI", "SppGEPE:Year.22.CI", "SppCHPE:Year.23.CI", "SppGEPE:Year.23.CI")

uncondCIMatrix<- matrix(NA, nrow=1, ncol=length(par.names), dimnames=list(c(1:length(1)), c(par.names.ci)))
uncondCIMatrix

uncondCIMatrix[1]<- uncondSEMatrix.filled[1]*1.96
uncondCIMatrix[2]<- uncondSEMatrix.filled[2]*1.96
uncondCIMatrix[3]<- uncondSEMatrix.filled[3]*1.96
uncondCIMatrix[4]<- uncondSEMatrix.filled[4]*1.96
uncondCIMatrix[5]<- uncondSEMatrix.filled[5]*1.96
uncondCIMatrix[6]<- uncondSEMatrix.filled[6]*1.96
uncondCIMatrix[7]<- uncondSEMatrix.filled[7]*1.96
uncondCIMatrix[8]<- uncondSEMatrix.filled[8]*1.96
uncondCIMatrix[9]<- uncondSEMatrix.filled[9]*1.96

uncondCIMatrix.filled<- uncondCIMatrix
uncondCIMatrix

t(uncondCIMatrix)

#-----
# Summary

AICtable

AICtable.2<- print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),digits = 5, LL = TRUE)

CoefMatrix

SEMatrix

t(ParalikMatrix.filled)

t(s.W.ParaEstMatrix.filled)

t(uncondSEMatrix.filled)

t(uncondCIMatrix)

#-----
# Variation in C/N between spp at PAL analyses. 5wk chicks, IDA analysis 4.

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see what files are available.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the data to be used.
DataFileName<- "Cks.5wk.PAL.INS.DataSet.Sub7.csv"

# Identify the models to be run.
ModelFileName<-"Cks.5wk.IsoNicheSpace.PAL.CandSet.C.csv" # 5 wk chk, PAL, C
ModelFileName<-"Cks.5wk.IsoNicheSpace.PAL.CandSet.N.csv" # 5 wk chk, PAL, N

#-----
# Import the data. Be sure to check the structure of the data and that R is reading continuous and categorical variables appropriately.
DataSet<-read.csv(DataFileName)

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

# Change parameter structures if needed.
DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

# Import candidate model set and convert to a character vector.
ModelSet.df<-read.csv(ModelFileName)  
ModelSet<-as.character(ModelSet.df$Model)

# View candidate models.
ModelSet

#-----
# Calculate all combinations of predictor variables to check candidate model set. Note, any candidate set may not include all combos of variables, but this will allow for checking that all combos have at least been considered (i.e., didnt forget some combo).

globalmodel<- lm(delta.15.N~Spp*Year.2*Spp:Year.2, data=DataSet)

AllCombCoefName<-names(globalmodel$coef)

# View all combinations of predictor variables.
AllCombCoefName

#-----
# Specify names of main explanatory parameters and associated SEs in the ModelSet (including interactions). Need to use : for interactions. !!Make sure that when creating any subsequent versions of CoefName in the code below that it follows the same order as is listed here. Further, when a candidate set has interaction models where (a) both main parameters are not included in the model, or where (b) the main effect is a categorical variable, you must include columns labeled here for all appropriate output!!
CoefName<- c("Intercept", "Intercept.SE", "SppCHPE", "SppCHPE.SE", "SppGEPE", "SppGEPE.SE", "Year.22", "Year.22.SE", "Year.23", "Year.23.SE", "SppCHPE:Year.22", "SppCHPE:Year.22.SE", "SppCHPE:Year.23", "SppCHPE:Year.23.SE", "SppGEPE:Year.22", "SppGEPE:Year.22.SE", "SppGEPE:Year.23", "SppGEPE:Year.23.SE")

# View predictor variables and associated SEs.
CoefName

#-----
# AIC function produces calculations for AIC Model Matrix that is used in conjunction with the ModelSet.
calculate.AIC<- function(AIC.Table, ModelSet) {
  
  deltaAIC.c<- AIC.Table$AIC.c-min(AIC.Table$AIC.c)
  lik.dAIC.c<- exp(-deltaAIC.c/2)
  AIC.c.W<- lik.dAIC.c/(sum(lik.dAIC.c))
  AICFinalMatrix<- data.frame(AICModelMatrix, deltaAIC.c, lik.dAIC.c, AIC.c.W)
}

#-----
# Create matrix to hold output from models.

# I. AIC Model Matrix, ncol=29 includes models, all coef names listed above, and column headings specified in the matrix below.
AICModelMatrix<- as.data.frame(matrix(NA, nrow=length(ModelSet), ncol=length(CoefName)+11, dimnames=list(c(1:length(ModelSet)), c("Models", CoefName, "N.Obs", "k", "EDF", "RMSE", "SSE", "logLik", "-2logLik", "mul.r.squared", "AIC", "AIC.c"))))

# View AIC Model Matrix.
head(AICModelMatrix)

#-----
# Loop for calculating model output and filling AIC Model Matrix. According to Burnham and Anderson, for least squares model fitting, K = total number of estimated regression parameters, including the intercept, and residual variation. In the case of R output, coef(model) includes an estimate of the intercept, thus, length(m$coef)+1 = K.

ModelOutput<- list()

for(i in 1:length(ModelSet)){
  
  ModelOutput[[i]]<- lm(as.formula(ModelSet[i]), data=DataSet)
  m<- ModelOutput[[i]]
  N.Obs<- nrow(DataSet)
  k<- length(m$coef)+1
  EDF<- N.Obs-k
  AIC<- AIC(m)
  AIC.c<- AIC+(2*k*(k+1))/(N.Obs-k-1)
  AICModelMatrix[i,"Models"]<- ModelSet[i]
  AICModelMatrix[i,"Intercept"]<- coef(m)["(Intercept)"]
  AICModelMatrix[i,"Intercept.SE"]<- summary(m)$coef["(Intercept)",2]
  AICModelMatrix[i,"SppCHPE"]<- coef(m)["SppCHPE"]
  AICModelMatrix[i,"SppCHPE.SE"]<- summary(m)$coef[,2]["SppCHPE"]
  AICModelMatrix[i,"SppGEPE"]<- coef(m)["SppGEPE"]
  AICModelMatrix[i,"SppGEPE.SE"]<- summary(m)$coef[,2]["SppGEPE"]
  AICModelMatrix[i,"Year.22"]<- coef(m)["Year.22"]
  AICModelMatrix[i,"Year.22.SE"]<- summary(m)$coef[,2]["Year.22"]
  AICModelMatrix[i,"Year.23"]<- coef(m)["Year.23"]
  AICModelMatrix[i,"Year.23.SE"]<- summary(m)$coef[,2]["Year.23"]
  AICModelMatrix[i,"SppCHPE:Year.22"]<- coef(m)["SppCHPE:Year.22"]
  AICModelMatrix[i,"SppCHPE:Year.22.SE"]<- summary(m)$coef[,2]["SppCHPE:Year.22"]
  AICModelMatrix[i,"SppCHPE:Year.23"]<- coef(m)["SppCHPE:Year.23"]
  AICModelMatrix[i,"SppCHPE:Year.23.SE"]<- summary(m)$coef[,2]["SppCHPE:Year.23"]
  AICModelMatrix[i,"SppGEPE:Year.22"]<- coef(m)["SppGEPE:Year.22"]
  AICModelMatrix[i,"SppGEPE:Year.22.SE"]<- summary(m)$coef[,2]["SppGEPE:Year.22"]
  AICModelMatrix[i,"SppGEPE:Year.23"]<- coef(m)["SppGEPE:Year.23"]
  AICModelMatrix[i,"SppGEPE:Year.23.SE"]<- summary(m)$coef[,2]["SppGEPE:Year.23"]
  AICModelMatrix[i,"N.Obs"]<- N.Obs
  AICModelMatrix[i,"k"]<- k 
  AICModelMatrix[i,"EDF"]<- EDF
  AICModelMatrix[i,"RMSE"]<- summary(m)$sigma
  AICModelMatrix[i,"SSE"]<- anova(m)["Residuals","Sum Sq"]
  AICModelMatrix[i,"logLik"]<- logLik(m)
  AICModelMatrix[i,"-2logLik"]<- -2*logLik(m)
  AICModelMatrix[i,"mul.r.squared"]<- summary(m)$r.squared
  AICModelMatrix[i,"AIC"]<- AIC
  AICModelMatrix[i,"AIC.c"]<- AIC.c
}

#-----
# Calculate deltaAICc, likdAICc, and AIC weights.
AIC.Output<-calculate.AIC(AICModelMatrix,as.character(ModelSet))

#-----
# View AIC Model Matrix.
print(AIC.Output)

# Write AIC Final Matrix output to .csv files.
write.table(AIC.Output, file="AICFinalMatrix.csv", col.names=NA, sep=",")

# Write model output to .csv files.
sink(paste("Model Summaries_",out="",".doc",sep=""))

for (i in 1:length(ModelSet)) {
  
  print(paste("MODEL ", i, sep = ""))
  print(ModelSet[[i]])
  print(summary(ModelOutput[[i]]))
  print(paste("-------------------------------------------------------------------"))
}

sink()

#-----
# Begin weighted and summed AIC calcs. Save previous file with AICFinalMatrix output, 
# excluding column A (Model #) and N.Obs through the lik.dAIC.c columns, 
# as weighted AIC Model Matrix Worksheet (W.AICFinalMatrixWS.csv). 

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see whats available. Make sure that new file created above is listed.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the new data WS to be used.
DataFileName2<- "W.AICFinalMatrixWS.csv"

# Import data.
DataSet2<-read.csv(DataFileName2)

# View DataSet. Note that all columns, with the exception of Models, should be listed as num.
head(DataSet2)
str(DataSet2)

#-----
# Specify a vector for DataSet2 Parameters. These will be used to define the parameters that will used to calculate Parameter likelihoods and W.ParaEsts that are multiplied by the model W. !!For interactions, be sure to specify them with . and not : as this is out they will be uploaded in DataSet2!!

# First recall CoefName and just exclude the .SE names in CoefName2 below.
CoefName<- c("Intercept", "Intercept.SE", "SppCHPE", "SppCHPE.SE", "SppGEPE", "SppGEPE.SE", "Year.22", "Year.22.SE", "Year.23", "Year.23.SE", "SppCHPE:Year.22", "SppCHPE:Year.22.SE", "SppCHPE:Year.23", "SppCHPE:Year.23.SE", "SppGEPE:Year.22", "SppGEPE:Year.22.SE", "SppGEPE:Year.23", "SppGEPE:Year.23.SE")

# View CoefName.
CoefName

CoefName2<- c("Intercept", "SppCHPE", "SppGEPE", "Year.22", "Year.23", "SppCHPE.Year.22", "SppCHPE.Year.23", "SppGEPE.Year.22", "SppGEPE.Year.23")

# View CoefName2.
CoefName2

#-----
# III. Calculate parameter likelihoods. Do this before W.ParaEst step below where DataSet2 NAs are turned into 0s.

# Create a vector to hold output.
Paralik<- vector()

# Create a matrix to hold output.
ParalikMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(CoefName2), dimnames=list(c(1:length(1)), c(CoefName2))))

# View ParalikMatrix.
ParalikMatrix

# Fill ParalikMatrix manually, b/c loop doesnt work with the subset function!

sub1<- subset(DataSet2, !is.na(Intercept))
nrow(sub1)
ParalikIntercept<- sum(sub1$AIC.c.W)
ParalikMatrix[1]<- signif(ParalikIntercept, digits=12)

sub2<- subset(DataSet2, !is.na(SppCHPE))
nrow(sub2)
ParalikSppCHPE<- sum(sub2$AIC.c.W)
ParalikMatrix[2]<- signif(ParalikSppCHPE, digits=12)

sub3<- subset(DataSet2, !is.na(SppGEPE))
nrow(sub3)
ParalikSppGEPE<- sum(sub3$AIC.c.W)
ParalikMatrix[3]<- signif(ParalikSppGEPE, digits=12)

sub4<- subset(DataSet2, !is.na(Year.22))
nrow(sub4)
ParalikYear.22<- sum(sub4$AIC.c.W)
ParalikMatrix[4]<- signif(ParalikYear.22, digits=12)

sub5<- subset(DataSet2, !is.na(Year.23))
nrow(sub5)
ParalikYear.23<- sum(sub5$AIC.c.W)
ParalikMatrix[5]<- signif(ParalikYear.23, digits=12)

sub6<- subset(DataSet2, !is.na(SppCHPE.Year.22))
nrow(sub6)
ParalikSppCHPE.Year.22<- sum(sub6$AIC.c.W)
ParalikMatrix[6]<- signif(ParalikSppCHPE.Year.22, digits=12)

sub7<- subset(DataSet2, !is.na(SppCHPE.Year.23))
nrow(sub7)
ParalikSppCHPE.Year.23<- sum(sub7$AIC.c.W)
ParalikMatrix[7]<- signif(ParalikSppCHPE.Year.23, digits=12)

sub8<- subset(DataSet2, !is.na(SppGEPE.Year.22))
nrow(sub8)
ParalikSppGEPE.Year.22<- sum(sub8$AIC.c.W)
ParalikMatrix[8]<- signif(ParalikSppGEPE.Year.22, digits=12)

sub9<- subset(DataSet2, !is.na(SppGEPE.Year.23))
nrow(sub9)
ParalikSppGEPE.Year.23<- sum(sub9$AIC.c.W)
ParalikMatrix[9]<- signif(ParalikSppGEPE.Year.23, digits=12)


#-----
# View ParalikMatrix
ParalikMatrix

#-----
# Transpose ParalikMatrix
transpose.ParalikMatrix<- t(ParalikMatrix)

# View transposed ParalikMatrix
transpose.ParalikMatrix

#-----
# Calculate Weighted Parameter Estimates (W.ParaEst).

# Turns NAs into 0s for calculating W.ParaEsts.
DataSet2[is.na(DataSet2)]<- 0

# View DataSet2 to check that NAs have have been replaced.
head(DataSet2)

# Specify a vector for new W.ParaEst.
W.CoefName2<- c("W.Intercept", "W.SppCHPE", "W.SppGEPE", "W.Year.22", "W.Year.23", "W.SppCHPE.Year.22", "W.SppCHPE.Year.23", "W.SppGEPE.Year.22", "W.SppGEPE.Year.23")

# View W.CoefName2.
W.CoefName2

# II. Create matrix to hold outputs for W.ParaEsts.

# Weighted ParaEst Matrix, ncol=8 includes all W.CoefName2 listed above.
W.ParaEstMatrix<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(W.CoefName2), dimnames=list(c(1:nrow(DataSet2)), c(W.CoefName2))))

# View W.ParaMatrix.
head(W.ParaEstMatrix)

# Run W.ParaEst loop to create W.ParaEsts by mulitplying the Para Ests for each model by the AIC.c weight of each model.

for(i in 1:length(CoefName2)) {
  
  W.ParaEst<- DataSet2[CoefName2[i]]*DataSet2["AIC.c.W"]
  W.ParaEstMatrix[,i]<- W.ParaEst
}

# View filled W.ParaEstMatrix.
W.ParaEstMatrix

#-----
# Summed W.ParaEsts.
# Specify a vector for summed W.ParaEsts.
s.W.CoefName2<- c("s.W.Intercept", "s.W.SppCHPE", "s.W.SppGEPE", "s.W.Year.22", "s.W.Year.23", "s.W.SppCHPE.Year.22", "s.W.SppCHPE.Year.23", "s.W.SppGEPE.Year.22", "s.W.SppGEPE.Year.23")

# View summed W.ParaEsts.
s.W.CoefName2

# Create 2 matrices to hold same output for summed W.ParaEsts

# II. 1-Summed W.ParaEsts Matrix, ncol=8 includes all s.W.CoefName2 listed above and nrow(DataSet2) to be used in the W.ParaEst.SE calcs.
s.W.ParaEstMatrix1<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(s.W.CoefName2), dimnames=list(c(1:nrow(DataSet2)), c(s.W.CoefName2))))

# 1-View s.W.ParaEstMatrix
s.W.ParaEstMatrix1

# 1-Run s.W.ParaEst loop. Calculates summed W.ParaEsts by summing all W.ParaEsts.

for(i in 1:length(W.CoefName2)) {
  
  s.W.ParaEst<- sum(W.ParaEstMatrix[W.CoefName2[i]])
  s.W.ParaEst[i]<- s.W.ParaEst
  s.W.ParaEstMatrix1[,i]<- s.W.ParaEst[i]
}

# View filled s.W.ParaEstMatrix1
s.W.ParaEstMatrix1

# III. 2-Summed W.ParaEsts Matrix, ncol=8 includes all s.W.CoefName2 listed above and nrow(1) to be used in final binding of all summed matrices.
s.W.ParaEstMatrix2<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName2), dimnames=list(c(1:length(1)), c(s.W.CoefName2))))

# 2-View s.W.ParaEstMatrix
s.W.ParaEstMatrix2

# 2-Run s.W.ParaEst loop. Calculates summed W.ParaEsts by summing all

for(i in 1:length(W.CoefName2)) {
  
  s.W.ParaEst<- sum(W.ParaEstMatrix[W.CoefName2[i]])
  s.W.ParaEst[i]<- s.W.ParaEst
  s.W.ParaEstMatrix2[,i]<- s.W.ParaEst[i]
}

# 2-View filled s.W.ParaEstMatrix2
s.W.ParaEstMatrix2

# 2-Transpose s.W.ParaEstMatrix2
transpose.s.W.ParaEstMatrix2<- t(s.W.ParaEstMatrix2)

# 2-View transposed s.W.ParaEstMatrix2
transpose.s.W.ParaEstMatrix2

#-----
# Calculate Weighted Parameter Estimate SEs (W.ParaEst.SE).

# Specify a vector for DataSet2 Parameter SEs. These will be used to define the parameter SEs that will used to calculate W.ParaEst.SE that are caculated below in the loop. !!Don't forget to use . instead of : for interactions as this is how they are uploaded in DataSet2.
CoefName3<- c("Intercept.SE", "SppCHPE.SE", "SppGEPE.SE", "Year.22.SE", "Year.23.SE", "SppCHPE.Year.22.SE", "SppCHPE.Year.23.SE", "SppGEPE.Year.22.SE", "SppGEPE.Year.23.SE")

# View CoefName3
CoefName3

# Specify a vector for new W.ParaEst.SE.
W.CoefName3.SE<- c("W.Intercept.SE", "W.SppCHPE.SE", "W.SppGEPE.SE", "W.Year.22.SE", "W.Year.23.SE", "W.SppCHPE.Year.22.SE", "W.SppCHPE.Year.23.SE", "W.SppGEPE.Year.22.SE", "W.SppGEPE.Year.23.SE")

# View W.CoefName2.SE
W.CoefName3.SE

# Create matrix to hold outputs for W.ParaEst.SE.

# II. Weighted ParaEst.SE Matrix, ncol=8 includes all W.CoefName3.SE listed above.
W.ParaEst.SEMatrix<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(W.CoefName3.SE), dimnames=list(c(1:nrow(DataSet2)), c(W.CoefName3.SE))))

# View W.ParaMatrix.
head(W.ParaEst.SEMatrix)

# Run W.ParaEst.SE loop to create W.ParaEstSEs using the following equation W.SE = (AIC.c.W*(sqrt((ParaEst SE^2)+(ParaEst-summedWParaEst)^2)))

for(j in 1:nrow(DataSet2)) {
  for(i in 1:length(CoefName3)) {
    
    if (DataSet2[j,CoefName3[i]] !=0) {W.ParaEst.SEMatrix[j,i]<- DataSet2$AIC.c.W[j]*(sqrt((DataSet2[j,CoefName3[i]])^2+(DataSet2[j,CoefName2[i]]-s.W.ParaEstMatrix1[j,i])^2))} else { 
      W.ParaEst.SEMatrix[j,i]<- 0}
  }
}

# View filled W.ParaEst.SEMatrix.
W.ParaEst.SEMatrix

#-----
# Summed W.ParaEst.SEs (Unconditional SEs)
# Specify a vector for summed W.ParaEsts.SEs.
s.W.CoefName3.SE<- c("s.W.Intercept.SE", "s.W.SppCHPE.SE", "s.W.SppGEPE.SE", "s.W.Year.22.SE", "s.W.Year.23.SE", "s.W.SppCHPE.Year.22.SE", "s.W.SppCHPE.Year.23.SE", "s.W.SppGEPE.Year.22.SE", "s.W.SppGEPE.Year.23.SE")

# View summed W.ParaEst.SEs.
s.W.CoefName3.SE

# Create matrix to hold outputs for summed W.ParaEst.SEs

# III. Summed W.ParaEst.SE Matrix, ncol=8 includes all s.W.CoefNames2.SE listed above
s.W.ParaEst.SEMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName3.SE), dimnames=list(c(1:length(1)), c(s.W.CoefName3.SE))))

# View s.W.ParaEst.SEMatrix
s.W.ParaEst.SEMatrix

#-----
# Run s.W.ParaEst.SE loop. Calculates summed W.ParaEsts by summing all

for(i in 1:length(s.W.CoefName3.SE)) {
  
  s.W.ParaEst.SE<- sum(W.ParaEst.SEMatrix[W.CoefName3.SE[i]])
  s.W.ParaEst.SE[i]<- s.W.ParaEst.SE
  s.W.ParaEst.SEMatrix[,i]<- s.W.ParaEst.SE[i]
}

# View filled s.W.ParaEst.SEMatrix
s.W.ParaEst.SEMatrix

# Transpose s.W.ParaEst.SEMatrix
transpose.s.W.ParaEst.SEMatrix<- t(s.W.ParaEst.SEMatrix)

# View transposed s.W.ParaEst.SEMatrix
transpose.s.W.ParaEst.SEMatrix

#-----
# Unconditional CIs
# Specify a vector for summed W.ParaEst.CIs.
s.W.CoefName3.CI<- c("s.W.Intercept.CI", "s.W.SppCHPE.CI", "s.W.SppGEPE.CI", "s.W.Year.22.CI", "s.W.Year.23.CI", "s.W.SppCHPE.Year.22.CI", "s.W.SppCHPE.Year.23.CI", "s.W.SppGEPE.Year.22.CI", "s.W.SppGEPE.Year.23.CI")

# View summed W.ParaEst.CIs.
s.W.CoefName3.CI

# Create matrix to hold outputs for summed W.ParaEst.CIs

# III. Summed W.ParaEst.CI Matrix, ncol=8 includes all s.W.CoefNames2.CI listed above
s.W.ParaEst.CIMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName3.CI), dimnames=list(c(1:length(1)), c(s.W.CoefName3.CI))))

# View s.W.ParaEst.CIMatrix
s.W.ParaEst.CIMatrix

# Run s.W.ParaEst.CI loop. Calculates summed W.ParaEst.CIs by multiplying 1.96*summed W.ParaEst.Ses

for(i in 1:length(s.W.CoefName3.SE)) {
  
  s.W.ParaEst.CI<- s.W.ParaEst.SEMatrix[s.W.CoefName3.SE[i]]*1.96
  s.W.ParaEst.CIMatrix[,i]<- s.W.ParaEst.CI
}

# View filled s.W.ParaEst.CIMatrix
s.W.ParaEst.CIMatrix

# Transpose s.W.ParaEst.CIMatrix
transpose.s.W.ParaEst.CIMatrix<- t(s.W.ParaEst.CIMatrix)

# View transposed s.W.ParaEst.CIMatrix
transpose.s.W.ParaEst.CIMatrix

#-----
# II. Bind DataSet2, W.ParaEstMatrix, W.ParaEst.SEMatrix

W.AICFinalMatrix<- cbind(DataSet2,W.ParaEstMatrix,W.ParaEst.SEMatrix)

# View W.AICFinaMatrix
W.AICFinalMatrix

# Write W.AICFinaMatrix to a .cvs file
write.table(W.AICFinalMatrix, file="W.AICFinalMatrix.csv", col.names=NA, sep=",")

#-----
# III. Bind ParalikMatrix, s.W.ParaEst.SEMatrix, s.W.ParaEst.CIMatrix
s.W.AICFinalMatrix<- cbind(transpose.ParalikMatrix,transpose.s.W.ParaEstMatrix2,transpose.s.W.ParaEst.SEMatrix,transpose.s.W.ParaEst.CIMatrix)

# View s.W.AICFinalMatrix
s.W.AICFinalMatrix

# Create vector for s.W.AICFinalMatrix column names.
s.W.AICFinalMatrixColName<- c("Paralik","s.W.ParaEst","UncondSE","UncondCI")

s.W.AICFinalMatrixColName

# Write s.W.AICFinaMatrix to a .cvs file. Note will have to move Column names over 1 column to be correct because R is weird.
write.table(s.W.AICFinalMatrix, file="s.W.AICFinalMatrix.csv", col.names=(s.W.AICFinalMatrixColName), sep=",")

#-----
# End code
# You should have now produced the following files with this code:
# 1-AICFinalMatrix.csv; use this file to sort models by deltaAIC.c and in excel calculate manually the cumAIC.c.W values.
# 2-Model Summaries_.doc; this is a word doc that holds the output from all models in the Candidate Set.
# 3-W.AICFinalMatrixWS.csv; this file is used as DataSet2 to run weighted calcs.
# 4-W.AICFinalMatrix.csv; this file contains all Models, Para Ests and associated SEs, AIC.c.Ws for each model, calculated W.Para Ests and associated W.SEs. These are the data that are used to create summed weighted estimates.
# 5-s.W.AICFinalMatrix.csv; this file holds summed weighted estimates including Parameter Likelihoods, s.W.ParaEsts, unconditional SEs, and unconditional CIs.
# !When all these files have been created, cut and paste them into 1 excel file with 4 sheets for each part of the analysis and label this excel file the same name as that for the Model Summaries_.doc file so that things can be referenced easily within the same analysis!


#-----
# Residual plots of most parameterized models.

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see what files are available.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the data to be used.
DataFileName<- "Cks.5wk.PAL.INS.DataSet.Sub6.csv"

#-----
# Import the data. Be sure to check the structure of the data and that R is reading continuous and categorical variables appropriately.
DataSet<-read.csv(DataFileName)

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

# Change parameter structures if needed.
DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

#-----
mcheck<- function(obj,...) {
  rs<- obj$resid
  fv<- obj$fitted
  par(mfrow=c(1,2))
  plot(fv,rs,xlab="Fitted Values",ylab="Residuals")
  abline(h=0,lty=2)
  qqnorm(rs,xlab="Normal scores",ylab="Ordered residuals",main="")
  par(mfrow=c(1,1))
  invisible(NULL)
}

#-----
# Residual plots of most parameterized model. Carbon.

DataSet

mod4<- lm(delta.13.C~Spp + Year.2 + Spp:Year.2, data=DataSet)
summary(mod4)

shapiro.test(resid(mod4))

#-----
# jpgs

jpeg(file="ResidPlot delta.13.C~Spp + Year.2 + Spp:Year.2.jpg")

mcheck(mod4)

dev.off()

jpeg(file="HistPlot.Resid delta.13.C~Spp + Year.2 + Spp:Year.2.jpg")

hist(resid(mod4))

dev.off()

#-----
# Nitrogen

mod1<- lm(delta.15.N~Spp + Year.2 + Spp:Year.2, data=DataSet)
summary(mod1)

shapiro.test(resid(mod1))

#-----
# jpgs

jpeg(file="ResidPlot delta.15.N~Spp + Year.2 + Spp:Year.2.jpg")

mcheck(mod1)

dev.off()

jpeg(file="HistPlot.Resid delta.15.N~Spp + Year.2 + Spp:Year.2.jpg")

hist(resid(mod1))

dev.off()

#-----
# Variation in C/N between spp across WAP analyses. 5wk chicks, IDA analysis 5.

# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see what files are available.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the data to be used.
DataFileName<- "Cks.5wk.WAP.INS.DataSet.Sub4.csv" #C
DataFileName<- "Cks.5wk.WAP.INS.DataSet.Sub3.csv" #N

# Identify the models to be run.
ModelFileName<-"Cks.5wk.IsoNicheSpace.WAP.CandSet.C.csv" #5wk Chk, WAP, C
ModelFileName<-"Cks.5wk.IsoNicheSpace.WAP.CandSet.N.csv" #5wk Chk, WAP, N

#-----
# Import the data. Be sure to check the structure of the data and that R is reading continuous and categorical variables appropriately.
DataSet<-read.csv(DataFileName)

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

# Change parameter structures if needed.
DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

# Import candidate model set and convert to a character vector.
ModelSet.df<-read.csv(ModelFileName)  
ModelSet<-as.character(ModelSet.df$Model)

# View candidate models.
ModelSet

#-----
# Calculate all combinations of predictor variables to check candidate model set. Note, any candidate set may not include all combos of variables, but this will allow for checking that all combos have at least been considered (i.e., didnt forget some combo).

globalmodel<- lm(delta.15.N~Location*Year.2*Location:Year.2, data=DataSet)

AllCombCoefName<-names(globalmodel$coef)

# View all combinations of predictor variables.
AllCombCoefName

#-----
# Specify names of main explanatory parameters and associated SEs in the ModelSet (including interactions). Need to use : for interactions. !!Make sure that when creating any subsequent versions of CoefName in the code below that it follows the same order as is listed here. Further, when a candidate set has interaction models where (a) both main parameters are not included in the model, or where (b) the main effect is a categorical variable, you must include columns labeled here for all appropriate output!!
CoefName<- c("Intercept", "Intercept.SE", "LocationCHA", "LocationCHA.SE", "LocationPAL", "LocationPAL.SE", "Year.22", "Year.22.SE", "Year.23", "Year.23.SE", "LocationPAL:Year.22", "LocationPAL:Year.22.SE", "LocationPAL:Year.23", "LocationPAL:Year.23.SE")

# View predictor variables and associated SEs.
CoefName

#-----
# AIC function produces calculations for AIC Model Matrix that is used in conjunction with the ModelSet.
calculate.AIC<- function(AIC.Table, ModelSet) {
  
  deltaAIC.c<- AIC.Table$AIC.c-min(AIC.Table$AIC.c)
  lik.dAIC.c<- exp(-deltaAIC.c/2)
  AIC.c.W<- lik.dAIC.c/(sum(lik.dAIC.c))
  AICFinalMatrix<- data.frame(AICModelMatrix, deltaAIC.c, lik.dAIC.c, AIC.c.W)
}

#-----
# Create matrix to hold output from models.

# I. AIC Model Matrix, ncol=29 includes models, all coef names listed above, and column headings specified in the matrix below.
AICModelMatrix<- as.data.frame(matrix(NA, nrow=length(ModelSet), ncol=length(CoefName)+11, dimnames=list(c(1:length(ModelSet)), c("Models", CoefName, "N.Obs", "k", "EDF", "RMSE", "SSE", "logLik", "-2logLik", "mul.r.squared", "AIC", "AIC.c"))))

# View AIC Model Matrix.
head(AICModelMatrix)

#-----
# Loop for calculating model output and filling AIC Model Matrix. According to Burnham and Anderson, for least squares model fitting, K = total number of estimated regression parameters, including the intercept, and residual variation. In the case of R output, coef(model) includes an estimate of the intercept, thus, length(m$coef)+1 = K.

ModelOutput<- list()

for(i in 1:length(ModelSet)){
  
  ModelOutput[[i]]<- lm(as.formula(ModelSet[i]), data=DataSet)
  m<- ModelOutput[[i]]
  N.Obs<- nrow(DataSet)
  k<- length(m$coef)+1
  EDF<- N.Obs-k
  AIC<- AIC(m)
  AIC.c<- AIC+(2*k*(k+1))/(N.Obs-k-1)
  AICModelMatrix[i,"Models"]<- ModelSet[i]
  AICModelMatrix[i,"Intercept"]<- coef(m)["(Intercept)"]
  AICModelMatrix[i,"Intercept.SE"]<- summary(m)$coef["(Intercept)",2]
  AICModelMatrix[i,"LocationCHA"]<- coef(m)["LocationCHA"]
  AICModelMatrix[i,"LocationCHA.SE"]<- summary(m)$coef[,2]["LocationCHA"]
  AICModelMatrix[i,"LocationPAL"]<- coef(m)["LocationPAL"]
  AICModelMatrix[i,"LocationPAL.SE"]<- summary(m)$coef[,2]["LocationPAL"]
  AICModelMatrix[i,"Year.22"]<- coef(m)["Year.22"]
  AICModelMatrix[i,"Year.22.SE"]<- summary(m)$coef[,2]["Year.22"]
  AICModelMatrix[i,"Year.23"]<- coef(m)["Year.23"]
  AICModelMatrix[i,"Year.23.SE"]<- summary(m)$coef[,2]["Year.23"]
  AICModelMatrix[i,"LocationPAL:Year.22"]<- coef(m)["LocationPAL:Year.22"]
  AICModelMatrix[i,"LocationPAL:Year.22.SE"]<- summary(m)$coef[,2]["LocationPAL:Year.22"]
  AICModelMatrix[i,"LocationPAL:Year.23"]<- coef(m)["LocationPAL:Year.23"]
  AICModelMatrix[i,"LocationPAL:Year.23.SE"]<- summary(m)$coef[,2]["LocationPAL:Year.23"]
  AICModelMatrix[i,"N.Obs"]<- N.Obs
  AICModelMatrix[i,"k"]<- k 
  AICModelMatrix[i,"EDF"]<- EDF
  AICModelMatrix[i,"RMSE"]<- summary(m)$sigma
  AICModelMatrix[i,"SSE"]<- anova(m)["Residuals","Sum Sq"]
  AICModelMatrix[i,"logLik"]<- logLik(m)
  AICModelMatrix[i,"-2logLik"]<- -2*logLik(m)
  AICModelMatrix[i,"mul.r.squared"]<- summary(m)$r.squared
  AICModelMatrix[i,"AIC"]<- AIC
  AICModelMatrix[i,"AIC.c"]<- AIC.c
}

#-----
# Calculate deltaAICc, likdAICc, and AIC weights.
AIC.Output<-calculate.AIC(AICModelMatrix,as.character(ModelSet))

#-----
# View AIC Model Matrix.
print(AIC.Output)

# Write AIC Final Matrix output to .csv files.
write.table(AIC.Output, file="AICFinalMatrix.csv", col.names=NA, sep=",")

# Write model output to .csv files.
sink(paste("Model Summaries_",out="",".doc",sep=""))

for (i in 1:length(ModelSet)) {
  
  print(paste("MODEL ", i, sep = ""))
  print(ModelSet[[i]])
  print(summary(ModelOutput[[i]]))
  print(paste("-------------------------------------------------------------------"))
}

sink()

#-----
# Begin weighted and summed AIC calcs. Save previous file with AICFinalMatrix output, excluding column A (Model #) and N.Obs through the lik.dAIC.c columns, as weighted AIC Model Matrix Worksheet (W.AICFinalMatrixWS.csv). For example, the scaup dataset, be sure to cut and paste the parameter output for "RFGinit:Rlipid", "RFGinit:Rlipid.SE" and paste into the column for "Rlipid:RFGinit", "Rlipid:RFGinit.SE".

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see whats available. Make sure that new file created above is listed.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the new data WS to be used.
DataFileName2<- "W.AICFinalMatrixWS.csv"

# Import data.
DataSet2<-read.csv(DataFileName2)

# View DataSet. Note that all columns, with the exception of Models, should be listed as num.
head(DataSet2)
str(DataSet2)

#-----
# Specify a vector for DataSet2 Parameters. These will be used to define the parameters that will used to calculate Parameter likelihoods and W.ParaEsts that are multiplied by the model W. !!For interactions, be sure to specify them with . and not : as this is out they will be uploaded in DataSet2!!

# First recall CoefName and just exclude the .SE names in CoefName2 below.
CoefName<- c("Intercept", "Intercept.SE", "LocationCHA", "LocationCHA.SE", "LocationPAL", "LocationPAL.SE", "Year.22", "Year.22.SE", "Year.23", "Year.23.SE", "LocationPAL.Year.22", "LocationPAL.Year.22.SE", "LocationPAL.Year.23", "LocationPAL.Year.23.SE")

# View CoefName.
CoefName

CoefName2<- c("Intercept", "LocationCHA", "LocationPAL", "Year.22", "Year.23", "LocationPAL.Year.22", "LocationPAL.Year.23")

# View CoefName2.
CoefName2

#-----
# III. Calculate parameter likelihoods. Do this before W.ParaEst step below where DataSet2 NAs are turned into 0s.

# Create a vector to hold output.
Paralik<- vector()

# Create a matrix to hold output.
ParalikMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(CoefName2), dimnames=list(c(1:length(1)), c(CoefName2))))

# View ParalikMatrix.
ParalikMatrix

# Fill ParalikMatrix manually, b/c loop doesnt work with the subset function!

sub1<- subset(DataSet2, !is.na(Intercept))
nrow(sub1)
ParalikIntercept<- sum(sub1$AIC.c.W)
ParalikMatrix[1]<- signif(ParalikIntercept, digits=12)

sub2<- subset(DataSet2, !is.na(LocationCHA))
nrow(sub2)
ParalikLocationCHA<- sum(sub2$AIC.c.W)
ParalikMatrix[2]<- signif(ParalikLocationCHA, digits=12)

sub3<- subset(DataSet2, !is.na(LocationPAL))
nrow(sub3)
ParalikLocationPAL<- sum(sub3$AIC.c.W)
ParalikMatrix[3]<- signif(ParalikLocationPAL, digits=12)

sub4<- subset(DataSet2, !is.na(Year.22))
nrow(sub4)
ParalikYear.22<- sum(sub4$AIC.c.W)
ParalikMatrix[4]<- signif(ParalikYear.22, digits=12)

sub5<- subset(DataSet2, !is.na(Year.23))
nrow(sub5)
ParalikYear.23<- sum(sub5$AIC.c.W)
ParalikMatrix[5]<- signif(ParalikYear.23, digits=12)

sub6<- subset(DataSet2, !is.na(LocationPAL.Year.22))
nrow(sub6)
ParalikLocationPAL.Year.22<- sum(sub6$AIC.c.W)
ParalikMatrix[6]<- signif(ParalikLocationPAL.Year.22, digits=12)

sub7<- subset(DataSet2, !is.na(LocationPAL.Year.23))
nrow(sub7)
ParalikLocationPAL.Year.23<- sum(sub7$AIC.c.W)
ParalikMatrix[7]<- signif(ParalikLocationPAL.Year.23, digits=12)

#-----
# View ParalikMatrix
ParalikMatrix

#-----
# Transpose ParalikMatrix
transpose.ParalikMatrix<- t(ParalikMatrix)

# View transposed ParalikMatrix
transpose.ParalikMatrix

#-----
# Calculate Weighted Parameter Estimates (W.ParaEst).

# Turns NAs into 0s for calculating W.ParaEsts.
DataSet2[is.na(DataSet2)]<- 0

# View DataSet2 to check that NAs have have been replaced.
head(DataSet2)

# Specify a vector for new W.ParaEst.
W.CoefName2<- c("W.Intercept", "W.LocationCHA", "W.LocationPAL", "W.Year.22", "W.Year.23", "W.LocationPAL.Year.22", "W.LocationPAL.Year.23")

# View W.CoefName2.
W.CoefName2

# II. Create matrix to hold outputs for W.ParaEsts.

# Weighted ParaEst Matrix, ncol=8 includes all W.CoefName2 listed above.
W.ParaEstMatrix<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(W.CoefName2), dimnames=list(c(1:nrow(DataSet2)), c(W.CoefName2))))

# View W.ParaMatrix.
head(W.ParaEstMatrix)

# Run W.ParaEst loop to create W.ParaEsts by mulitplying the Para Ests for each model by the AIC.c weight of each model.

for(i in 1:length(CoefName2)) {
  
  W.ParaEst<- DataSet2[CoefName2[i]]*DataSet2["AIC.c.W"]
  W.ParaEstMatrix[,i]<- W.ParaEst
}

# View filled W.ParaEstMatrix.
W.ParaEstMatrix

#-----
# Summed W.ParaEsts.
# Specify a vector for summed W.ParaEsts.
s.W.CoefName2<- c("s.W.Intercept", "s.W.LocationCHA", "s.W.LocationPAL", "s.W.Year.22", "s.W.Year.23", "s.W.LocationPAL.Year.22", "s.W.LocationPAL.Year.23")

# View summed W.ParaEsts.
s.W.CoefName2

# Create 2 matrices to hold same output for summed W.ParaEsts

# II. 1-Summed W.ParaEsts Matrix, ncol=8 includes all s.W.CoefName2 listed above and nrow(DataSet2) to be used in the W.ParaEst.SE calcs.
s.W.ParaEstMatrix1<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(s.W.CoefName2), dimnames=list(c(1:nrow(DataSet2)), c(s.W.CoefName2))))

# 1-View s.W.ParaEstMatrix
s.W.ParaEstMatrix1

# 1-Run s.W.ParaEst loop. Calculates summed W.ParaEsts by summing all W.ParaEsts.

for(i in 1:length(W.CoefName2)) {
  
  s.W.ParaEst<- sum(W.ParaEstMatrix[W.CoefName2[i]])
  s.W.ParaEst[i]<- s.W.ParaEst
  s.W.ParaEstMatrix1[,i]<- s.W.ParaEst[i]
}

# View filled s.W.ParaEstMatrix1
s.W.ParaEstMatrix1

# III. 2-Summed W.ParaEsts Matrix, ncol=8 includes all s.W.CoefName2 listed above and nrow(1) to be used in final binding of all summed matrices.
s.W.ParaEstMatrix2<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName2), dimnames=list(c(1:length(1)), c(s.W.CoefName2))))

# 2-View s.W.ParaEstMatrix
s.W.ParaEstMatrix2

# 2-Run s.W.ParaEst loop. Calculates summed W.ParaEsts by summing all

for(i in 1:length(W.CoefName2)) {
  
  s.W.ParaEst<- sum(W.ParaEstMatrix[W.CoefName2[i]])
  s.W.ParaEst[i]<- s.W.ParaEst
  s.W.ParaEstMatrix2[,i]<- s.W.ParaEst[i]
}

# 2-View filled s.W.ParaEstMatrix2
s.W.ParaEstMatrix2

# 2-Transpose s.W.ParaEstMatrix2
transpose.s.W.ParaEstMatrix2<- t(s.W.ParaEstMatrix2)

# 2-View transposed s.W.ParaEstMatrix2
transpose.s.W.ParaEstMatrix2

#-----
# Calculate Weighted Parameter Estimate SEs (W.ParaEst.SE).

# Specify a vector for DataSet2 Parameter SEs. These will be used to define the parameter SEs that will used to calculate W.ParaEst.SE that are caculated below in the loop. !!Don't forget to use . instead of : for interactions as this is how they are uploaded in DataSet2.
CoefName3<- c("Intercept.SE", "LocationCHA.SE", "LocationPAL.SE", "Year.22.SE", "Year.23.SE", "LocationPAL.Year.22.SE", "LocationPAL.Year.23.SE")

# View CoefName3
CoefName3

# Specify a vector for new W.ParaEst.SE.
W.CoefName3.SE<- c("W.Intercept.SE", "W.LocationCHA.SE", "W.LocationPAL.SE", "W.Year.22.SE", "W.Year.23.SE", "W.LocationPAL.Year.22.SE", "W.LocationPAL.Year.23.SE")

# View W.CoefName2.SE
W.CoefName3.SE

# Create matrix to hold outputs for W.ParaEst.SE.

# II. Weighted ParaEst.SE Matrix, ncol=8 includes all W.CoefName3.SE listed above.
W.ParaEst.SEMatrix<- as.data.frame(matrix(NA, nrow=nrow(DataSet2), ncol=length(W.CoefName3.SE), dimnames=list(c(1:nrow(DataSet2)), c(W.CoefName3.SE))))

# View W.ParaMatrix.
head(W.ParaEst.SEMatrix)

# Run W.ParaEst.SE loop to create W.ParaEstSEs using the following equation W.SE = (AIC.c.W*(sqrt((ParaEst SE^2)+(ParaEst-summedWParaEst)^2)))

for(j in 1:nrow(DataSet2)) {
  for(i in 1:length(CoefName3)) {
    
    if (DataSet2[j,CoefName3[i]] !=0) {W.ParaEst.SEMatrix[j,i]<- DataSet2$AIC.c.W[j]*(sqrt((DataSet2[j,CoefName3[i]])^2+(DataSet2[j,CoefName2[i]]-s.W.ParaEstMatrix1[j,i])^2))} else { 
      W.ParaEst.SEMatrix[j,i]<- 0}
  }
}

# View filled W.ParaEst.SEMatrix.
W.ParaEst.SEMatrix

#-----
# Summed W.ParaEst.SEs (Unconditional SEs)
# Specify a vector for summed W.ParaEsts.SEs.
s.W.CoefName3.SE<- c("s.W.Intercept.SE", "s.W.LocationCHA.SE", "s.W.LocationPAL.SE", "s.W.Year.22.SE", "s.W.Year.23.SE", "s.W.LocationPAL.Year.22.SE", "s.W.LocationPAL.Year.23.SE")

# View summed W.ParaEst.SEs.
s.W.CoefName3.SE

# Create matrix to hold outputs for summed W.ParaEst.SEs

# III. Summed W.ParaEst.SE Matrix, ncol=8 includes all s.W.CoefNames2.SE listed above
s.W.ParaEst.SEMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName3.SE), dimnames=list(c(1:length(1)), c(s.W.CoefName3.SE))))

# View s.W.ParaEst.SEMatrix
s.W.ParaEst.SEMatrix

#-----
# Run s.W.ParaEst.SE loop. Calculates summed W.ParaEsts by summing all

for(i in 1:length(s.W.CoefName3.SE)) {
  
  s.W.ParaEst.SE<- sum(W.ParaEst.SEMatrix[W.CoefName3.SE[i]])
  s.W.ParaEst.SE[i]<- s.W.ParaEst.SE
  s.W.ParaEst.SEMatrix[,i]<- s.W.ParaEst.SE[i]
}

# View filled s.W.ParaEst.SEMatrix
s.W.ParaEst.SEMatrix

# Transpose s.W.ParaEst.SEMatrix
transpose.s.W.ParaEst.SEMatrix<- t(s.W.ParaEst.SEMatrix)

# View transposed s.W.ParaEst.SEMatrix
transpose.s.W.ParaEst.SEMatrix

#-----
# Unconditional CIs
# Specify a vector for summed W.ParaEst.CIs.
s.W.CoefName3.CI<- c("s.W.Intercept.CI", "s.W.LocationCHA.CI", "s.W.LocationPAL.CI", "s.W.Year.22.CI", "s.W.Year.23.CI", "s.W.LocationPAL.Year.22.CI", "s.W.LocationPAL.Year.23.CI")

# View summed W.ParaEst.CIs.
s.W.CoefName3.CI

# Create matrix to hold outputs for summed W.ParaEst.CIs

# III. Summed W.ParaEst.CI Matrix, ncol=8 includes all s.W.CoefNames2.CI listed above
s.W.ParaEst.CIMatrix<- as.data.frame(matrix(NA, nrow=1, ncol=length(s.W.CoefName3.CI), dimnames=list(c(1:length(1)), c(s.W.CoefName3.CI))))

# View s.W.ParaEst.CIMatrix
s.W.ParaEst.CIMatrix

# Run s.W.ParaEst.CI loop. Calculates summed W.ParaEst.CIs by multiplying 1.96*summed W.ParaEst.Ses

for(i in 1:length(s.W.CoefName3.SE)) {
  
  s.W.ParaEst.CI<- s.W.ParaEst.SEMatrix[s.W.CoefName3.SE[i]]*1.96
  s.W.ParaEst.CIMatrix[,i]<- s.W.ParaEst.CI
}

# View filled s.W.ParaEst.CIMatrix
s.W.ParaEst.CIMatrix

# Transpose s.W.ParaEst.CIMatrix
transpose.s.W.ParaEst.CIMatrix<- t(s.W.ParaEst.CIMatrix)

# View transposed s.W.ParaEst.CIMatrix
transpose.s.W.ParaEst.CIMatrix

#-----
# II. Bind DataSet2, W.ParaEstMatrix, W.ParaEst.SEMatrix

W.AICFinalMatrix<- cbind(DataSet2,W.ParaEstMatrix,W.ParaEst.SEMatrix)

# View W.AICFinaMatrix
W.AICFinalMatrix

# Write W.AICFinaMatrix to a .cvs file
write.table(W.AICFinalMatrix, file="W.AICFinalMatrix.csv", col.names=NA, sep=",")

#-----
# III. Bind ParalikMatrix, s.W.ParaEst.SEMatrix, s.W.ParaEst.CIMatrix
s.W.AICFinalMatrix<- cbind(transpose.ParalikMatrix,transpose.s.W.ParaEstMatrix2,transpose.s.W.ParaEst.SEMatrix,transpose.s.W.ParaEst.CIMatrix)

# View s.W.AICFinalMatrix
s.W.AICFinalMatrix

# Create vector for s.W.AICFinalMatrix column names.
s.W.AICFinalMatrixColName<- c("Paralik","s.W.ParaEst","UncondSE","UncondCI")

s.W.AICFinalMatrixColName

# Write s.W.AICFinaMatrix to a .cvs file. Note will have to move Column names over 1 column to be correct because R is weird.
write.table(s.W.AICFinalMatrix, file="s.W.AICFinalMatrix.csv", col.names=(s.W.AICFinalMatrixColName), sep=",")

#-----
# End code.
# You should have now produced the following files with this code:
# 1-AICFinalMatrix.csv; use this file to sort models by deltaAIC.c and in excel calculate manually the cumAIC.c.W values.
# 2-Model Summaries_.doc; this is a word doc that holds the output from all models in the Candidate Set.
# 3-W.AICFinalMatrixWS.csv; this file is used as DataSet2 to run weighted calcs.
# 4-W.AICFinalMatrix.csv; this file contains all Models, Para Ests and associated SEs, AIC.c.Ws for each model, calculated W.Para Ests and associated W.SEs. These are the data that are used to create summed weighted estimates.
# 5-s.W.AICFinalMatrix.csv; this file holds summed weighted estimates including Parameter Likelihoods, s.W.ParaEsts, unconditional SEs, and unconditional CIs.
# !When all these files have been created, cut and paste them into 1 excel file with 4 sheets for each part of the analysis and label this excel file the same name as that for the Model Summaries_.doc file so that things can be referenced easily within the same analysis!

#-----
# Residual plots of most parameterized models.

#-----
# Reset R, clear all objects.
rm(list=ls())

#-----
# Specify working directory and see what files are available.
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

# Identify the data to be used.
DataFileName<- "Cks.5wk.WAP.INS.DataSet.Sub4.csv" #C
DataFileName<- "Cks.5wk.WAP.INS.DataSet.Sub3.csv" #N

#-----
# Import the data. Be sure to check the structure of the data and that R is reading continuous and categorical variables appropriately.
DataSet<-read.csv(DataFileName)

# View DataSet.
head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

# Change parameter structures if needed.
DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet


#-----
mcheck<- function(obj,...) {
  rs<- obj$resid
  fv<- obj$fitted
  par(mfrow=c(1,2))
  plot(fv,rs,xlab="Fitted Values",ylab="Residuals")
  abline(h=0,lty=2)
  qqnorm(rs,xlab="Normal scores",ylab="Ordered residuals",main="")
  par(mfrow=c(1,1))
  invisible(NULL)
}

#-----
# Residual plots of most parameterized model. Carbon.

DataSet

mod4<- lm(delta.13.C~Location + Year.2 + Location:Year.2, data=DataSet)
summary(mod4)

shapiro.test(resid(mod4))

#-----
# jpgs

jpeg(file="ResidPlot delta.13.C~Location + Year.2 + Location:Year.2.jpg")

mcheck(mod4)

dev.off()

jpeg(file="HistPlot.Resid delta.13.C~Location + Year.2 + Location:Year.2.jpg")

hist(resid(mod4))

dev.off()

#-----
# Nitrogen

mod1<- lm(delta.15.N~Location + Year.2 + Location:Year.2, data=DataSet)
summary(mod1)

shapiro.test(resid(mod1))

#-----
# jpgs

#-----
jpeg(file="ResidPlot delta.15.N~Location + Year.2 + Location:Year.2.jpg")

mcheck(mod1)

dev.off()

jpeg(file="HistPlot.Resid delta.15.N~Location + Year.2 + Location:Year.2.jpg")

hist(resid(mod1))

dev.off()

#-----
#Figures 4-6
# Species isotope graphs by year. Adults.

#-----
# Reset R, clear all objects
rm(list=ls()) 

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

DataSet<- read.csv("Ads.INS.DataSet.Sub3_FINAL.csv")

head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

range(DataSet$delta.15.N)

#-----
DataSet1<- DataSet

head(DataSet1)
summary(DataSet1)
nrow(DataSet1)
str(DataSet1)

#-----
list.files()

DataSet2<-read.csv("Ads.INS.DataSet.Sub4_FINAL.csv")

head(DataSet2)
summary(DataSet2)
nrow(DataSet2)
str(DataSet2)

DataSet2$Year.2<- as.factor(DataSet2$Year.2)
str(DataSet2)
DataSet2

range(DataSet2$delta.13.C)

#-----
# Subset DataSets into spp and year, by isotope

#-----
DataSet1

N.Sub1<- subset(DataSet1,Year.2=="1")
N.Sub2<- subset(DataSet1,Year.2=="2")
N.Sub3<- subset(DataSet1,Year.2=="3") 

N.Sub1.Spp1<- subset(N.Sub1,Spp=="ADPE")
N.Sub1.Spp2<- subset(N.Sub1,Spp=="CHPE")
N.Sub1.Spp3<- subset(N.Sub1,Spp=="GEPE")

N.Sub2.Spp1<- subset(N.Sub2,Spp=="ADPE")
N.Sub2.Spp2<- subset(N.Sub2,Spp=="CHPE")
N.Sub2.Spp3<- subset(N.Sub2,Spp=="GEPE")

N.Sub3.Spp1<- subset(N.Sub3,Spp=="ADPE")
N.Sub3.Spp2<- subset(N.Sub3,Spp=="CHPE")
N.Sub3.Spp3<- subset(N.Sub3,Spp=="GEPE")

#-----
DataSet2

C.Sub1<- subset(DataSet2,Year.2=="1")
C.Sub2<- subset(DataSet2,Year.2=="2")
C.Sub3<- subset(DataSet2,Year.2=="3") 

C.Sub1.Spp1<- subset(C.Sub1,Spp=="ADPE")
C.Sub1.Spp2<- subset(C.Sub1,Spp=="CHPE")
C.Sub1.Spp3<- subset(C.Sub1,Spp=="GEPE")

C.Sub2.Spp1<- subset(C.Sub2,Spp=="ADPE")
C.Sub2.Spp2<- subset(C.Sub2,Spp=="CHPE")
C.Sub2.Spp3<- subset(C.Sub2,Spp=="GEPE")

C.Sub3.Spp1<- subset(C.Sub3,Spp=="ADPE")
C.Sub3.Spp2<- subset(C.Sub3,Spp=="CHPE")
C.Sub3.Spp3<- subset(C.Sub3,Spp=="GEPE")

#-----
# Averages for each spp by year for each isotope, and standard error. The variance estimate (var) is a built in function in R. The stdev is the sqrt of the variance. The se function was checked in excel and gave the same value.

variance<- function(x) sum((x-mean(x))^2)/(length(x)-1)
stdev<- function(x) sqrt(sum((x-mean(x))^2)/(length(x)-1))
se<- function(x) sqrt(var(x)/length(x))

# ADPE

AP.Yr1.mod1<- mean(C.Sub1.Spp1$delta.13.C)
AP.Yr1.mod1.se<- se(C.Sub1.Spp1$delta.13.C)

AP.Yr1.mod2<- mean(N.Sub1.Spp1$delta.15.N)
AP.Yr1.mod2.se<- se(N.Sub1.Spp1$delta.15.N)

AP.Yr2.mod3<- mean(C.Sub2.Spp1$delta.13.C)
AP.Yr2.mod3.se<- se(C.Sub2.Spp1$delta.13.C)

AP.Yr2.mod4<- mean(N.Sub2.Spp1$delta.15.N)
AP.Yr2.mod4.se<- se(N.Sub2.Spp1$delta.15.N)

AP.Yr3.mod5<- mean(C.Sub3.Spp1$delta.13.C)
AP.Yr3.mod5.se<- se(C.Sub3.Spp1$delta.13.C)

AP.Yr3.mod6<- mean(N.Sub3.Spp1$delta.15.N)
AP.Yr3.mod6.se<- se(N.Sub3.Spp1$delta.15.N)

#-----
# CHPE

CP.Yr1.mod1<- mean(C.Sub1.Spp2$delta.13.C)
CP.Yr1.mod1.se<- se(C.Sub1.Spp2$delta.13.C)

CP.Yr1.mod2<- mean(N.Sub1.Spp2$delta.15.N)
CP.Yr1.mod2.se<- se(N.Sub1.Spp2$delta.15.N)

CP.Yr2.mod3<- mean(C.Sub2.Spp2$delta.13.C)
CP.Yr2.mod3.se<- se(C.Sub2.Spp2$delta.13.C)

CP.Yr2.mod4<- mean(N.Sub2.Spp2$delta.15.N)
CP.Yr2.mod4.se<- se(N.Sub2.Spp2$delta.15.N)

CP.Yr3.mod5<- mean(C.Sub3.Spp2$delta.13.C)
CP.Yr3.mod5.se<- se(C.Sub3.Spp2$delta.13.C)

CP.Yr3.mod6<- mean(N.Sub3.Spp2$delta.15.N)
CP.Yr3.mod6.se<- se(N.Sub3.Spp2$delta.15.N)

#-----
# GEPE

GP.Yr1.mod1<- mean(C.Sub1.Spp3$delta.13.C)
GP.Yr1.mod1.se<- se(C.Sub1.Spp3$delta.13.C)

GP.Yr1.mod2<- mean(N.Sub1.Spp3$delta.15.N)
GP.Yr1.mod2.se<- se(N.Sub1.Spp3$delta.15.N)

GP.Yr2.mod3<- mean(C.Sub2.Spp3$delta.13.C)
GP.Yr2.mod3.se<- se(C.Sub2.Spp3$delta.13.C)

GP.Yr2.mod4<- mean(N.Sub2.Spp3$delta.15.N)
GP.Yr2.mod4.se<- se(N.Sub2.Spp3$delta.15.N)

GP.Yr3.mod5<- mean(C.Sub3.Spp3$delta.13.C)
GP.Yr3.mod5.se<- se(C.Sub3.Spp3$delta.13.C)

GP.Yr3.mod6<- mean(N.Sub3.Spp3$delta.15.N)
GP.Yr3.mod6.se<- se(N.Sub3.Spp3$delta.15.N)

#-----
# ADPE. Plot each year by isotopes

AP.x1<- AP.Yr1.mod1
AP.y1<- AP.Yr1.mod2
AP.x1.se<- AP.Yr1.mod1.se
AP.y1.se<- AP.Yr1.mod2.se

AP.x2<- AP.Yr2.mod3
AP.y2<- AP.Yr2.mod4
AP.x2.se<- AP.Yr2.mod3.se
AP.y2.se<- AP.Yr2.mod4.se

AP.x3<- AP.Yr3.mod5
AP.y3<- AP.Yr3.mod6
AP.x3.se<- AP.Yr3.mod5.se
AP.y3.se<- AP.Yr3.mod6.se

#-----
# CHPE. Plot each year by isotopes

CP.x1<- CP.Yr1.mod1
CP.y1<- CP.Yr1.mod2
CP.x1.se<- CP.Yr1.mod1.se
CP.y1.se<- CP.Yr1.mod2.se

CP.x2<- CP.Yr2.mod3
CP.y2<- CP.Yr2.mod4
CP.x2.se<- CP.Yr2.mod3.se
CP.y2.se<- CP.Yr2.mod4.se

CP.x3<- CP.Yr3.mod5
CP.y3<- CP.Yr3.mod6
CP.x3.se<- CP.Yr3.mod5.se
CP.y3.se<- CP.Yr3.mod6.se

#-----
# GEPE. Plot each year by isotopes

GP.x1<- GP.Yr1.mod1
GP.y1<- GP.Yr1.mod2
GP.x1.se<- GP.Yr1.mod1.se
GP.y1.se<- GP.Yr1.mod2.se

GP.x2<- GP.Yr2.mod3
GP.y2<- GP.Yr2.mod4
GP.x2.se<- GP.Yr2.mod3.se
GP.y2.se<- GP.Yr2.mod4.se

GP.x3<- GP.Yr3.mod5
GP.y3<- GP.Yr3.mod6
GP.x3.se<- GP.Yr3.mod5.se
GP.y3.se<- GP.Yr3.mod6.se

#-----
# Plot

tiff(file="ADs.CN.Year3.tiff",res=300,width=15,height=15,unit="cm")
par(mfrow=c(1,3))

plot1<- plot(AP.y1~AP.x1, main="2007/08", cex.main=1.5, pch=19, col="black", xlim=c(-28,-24.0), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y1~CP.x1, pch=19, col="orangered1")
points(GP.y1~GP.x1, pch=19, col="darkblue")

points(AP.x1, AP.y1, pch=19, col="black", bg="black")
arrows(AP.x1, AP.y1-AP.y1.se, AP.x1, AP.y1+AP.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x1-AP.x1.se, AP.y1, AP.x1+AP.x1.se, AP.y1, code=3, angle=90, length=0.1, col="black")

points(CP.x1, CP.y1, pch=19, col="orangered", bg="orangered")
arrows(CP.x1, CP.y1-CP.y1.se, CP.x1, CP.y1+CP.y1.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x1-CP.x1.se, CP.y1, CP.x1+CP.x1.se, CP.y1, code=3, angle=90, length=0.1, col="orangered")

points(GP.x1, GP.y1, pch=19, col="darkblue", bg="darkblue")
arrows(GP.x1, GP.y1-GP.y1.se, GP.x1, GP.y1+GP.y1.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x1-GP.x1.se, GP.y1, GP.x1+GP.x1.se, GP.y1, code=3, angle=90, length=0.1, col="darkblue")

plot2<- plot(AP.y2~AP.x2, main="2008/09", cex.main=1.5, pch=19, col="black", xlim=c(-28,-24.0), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y2~CP.x2, pch=19, col="orangered1")
points(GP.y2~GP.x2, pch=19, col="darkblue")

points(AP.x2, AP.y2, pch=19, col="black", bg="black")
arrows(AP.x2, AP.y2-AP.y2.se, AP.x2, AP.y2+AP.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x2-AP.x2.se, AP.y2, AP.x2+AP.x2.se, AP.y2, code=3, angle=90, length=0.1, col="black")

points(CP.x2, CP.y2, pch=19, col="orangered", bg="orangered")
arrows(CP.x2, CP.y2-CP.y2.se, CP.x2, CP.y2+CP.y2.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x2-CP.x2.se, CP.y2, CP.x2+CP.x2.se, CP.y2, code=3, angle=90, length=0.1, col="orangered")

points(GP.x2, GP.y2, pch=19, col="darkblue", bg="darkblue")
arrows(GP.x2, GP.y2-GP.y2.se, GP.x2, GP.y2+GP.y2.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x2-GP.x2.se, GP.y2, GP.x2+GP.x2.se, GP.y2, code=3, angle=90, length=0.1, col="darkblue")

plot3<- plot(AP.y3~AP.x3, main="2009/10", cex.main=1.5, pch=19, col="black", xlim=c(-28,-24.0), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y3~CP.x3, pch=19, col="orangered1")
points(GP.y3~GP.x3, pch=19, col="darkblue")

points(AP.x3, AP.y3, pch=19, col="black", bg="black")
arrows(AP.x3, AP.y3-AP.y3.se, AP.x3, AP.y3+AP.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x3-AP.x3.se, AP.y3, AP.x3+AP.x3.se, AP.y3, code=3, angle=90, length=0.1, col="black")

points(CP.x3, CP.y3, pch=19, col="orangered", bg="orangered")
arrows(CP.x3, CP.y3-CP.y3.se, CP.x3, CP.y3+CP.y3.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x3-CP.x3.se, CP.y3, CP.x3+CP.x3.se, CP.y3, code=3, angle=90, length=0.1, col="orangered")

points(GP.x3, GP.y3, pch=19, col="darkblue", bg="darkblue")
arrows(GP.x3, GP.y3-GP.y3.se, GP.x3, GP.y3+GP.y3.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x3-GP.x3.se, GP.y3, GP.x3+GP.x3.se, GP.y3, code=3, angle=90, length=0.1, col="darkblue")

dev.off()

#-----
# Species isotope graphs by year, Cks.D5.PAL

#-----
# Reset R, clear all objects
rm(list=ls()) 

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

DataSet<- read.csv("Cks.D5.PAL.INS.DataSet.Sub2_FINAL.csv")

head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

DataSet$Year.2<- as.factor(DataSet$Year.2)
DataSet$Nest<- as.factor(DataSet$Nest)
str(DataSet)
DataSet

range(DataSet$delta.15.N)
range(DataSet$delta.13.C)

#-----
DataSet1<- DataSet

head(DataSet1)
summary(DataSet1)
nrow(DataSet1)
str(DataSet1)

#-----
# Subset DataSets into spp and year, by isotope. Can use same dataset for both N and C.

#-----
# DataSet1.N

N.Sub1<- subset(DataSet1,Year.2=="1")
N.Sub2<- subset(DataSet1,Year.2=="2")
N.Sub3<- subset(DataSet1,Year.2=="3") 

N.Sub1.Spp1<- subset(N.Sub1,Spp=="ADPE")
N.Sub1.Spp2<- subset(N.Sub1,Spp=="CHPE")
N.Sub1.Spp3<- subset(N.Sub1,Spp=="GEPE")

N.Sub2.Spp1<- subset(N.Sub2,Spp=="ADPE")
N.Sub2.Spp2<- subset(N.Sub2,Spp=="CHPE")
N.Sub2.Spp3<- subset(N.Sub2,Spp=="GEPE")

N.Sub3.Spp1<- subset(N.Sub3,Spp=="ADPE")
N.Sub3.Spp2<- subset(N.Sub3,Spp=="CHPE")
N.Sub3.Spp3<- subset(N.Sub3,Spp=="GEPE")

#-----
# DataSet1.C

C.Sub1<- subset(DataSet1,Year.2=="1")
C.Sub2<- subset(DataSet1,Year.2=="2")
C.Sub3<- subset(DataSet1,Year.2=="3") 

C.Sub1.Spp1<- subset(C.Sub1,Spp=="ADPE")
C.Sub1.Spp2<- subset(C.Sub1,Spp=="CHPE")
C.Sub1.Spp3<- subset(C.Sub1,Spp=="GEPE")

C.Sub2.Spp1<- subset(C.Sub2,Spp=="ADPE")
C.Sub2.Spp2<- subset(C.Sub2,Spp=="CHPE")
C.Sub2.Spp3<- subset(C.Sub2,Spp=="GEPE")

C.Sub3.Spp1<- subset(C.Sub3,Spp=="ADPE")
C.Sub3.Spp2<- subset(C.Sub3,Spp=="CHPE")
C.Sub3.Spp3<- subset(C.Sub3,Spp=="GEPE")

#-----
# Averages for each spp by year for each isotope, and standard error. The variance estimate (var) is a built in function in R. The stdev is the sqrt of the variance. The se function was checked in excel and gave the same value.

variance<- function(x) sum((x-mean(x))^2)/(length(x)-1)
stdev<- function(x) sqrt(sum((x-mean(x))^2)/(length(x)-1))
se<- function(x) sqrt(var(x)/length(x))

# ADPE

AP.Yr1.mod1<- mean(C.Sub1.Spp1$delta.13.C)
AP.Yr1.mod1.se<- se(C.Sub1.Spp1$delta.13.C)

AP.Yr1.mod2<- mean(N.Sub1.Spp1$delta.15.N)
AP.Yr1.mod2.se<- se(N.Sub1.Spp1$delta.15.N)

AP.Yr2.mod3<- mean(C.Sub2.Spp1$delta.13.C)
AP.Yr2.mod3.se<- se(C.Sub2.Spp1$delta.13.C)

AP.Yr2.mod4<- mean(N.Sub2.Spp1$delta.15.N)
AP.Yr2.mod4.se<- se(N.Sub2.Spp1$delta.15.N)

AP.Yr3.mod5<- mean(C.Sub3.Spp1$delta.13.C)
AP.Yr3.mod5.se<- se(C.Sub3.Spp1$delta.13.C)

AP.Yr3.mod6<- mean(N.Sub3.Spp1$delta.15.N)
AP.Yr3.mod6.se<- se(N.Sub3.Spp1$delta.15.N)

#-----
# CHPE

CP.Yr1.mod1<- mean(C.Sub1.Spp2$delta.13.C)
CP.Yr1.mod1.se<- se(C.Sub1.Spp2$delta.13.C)

CP.Yr1.mod2<- mean(N.Sub1.Spp2$delta.15.N)
CP.Yr1.mod2.se<- se(N.Sub1.Spp2$delta.15.N)

CP.Yr2.mod3<- mean(C.Sub2.Spp2$delta.13.C)
CP.Yr2.mod3.se<- se(C.Sub2.Spp2$delta.13.C)

CP.Yr2.mod4<- mean(N.Sub2.Spp2$delta.15.N)
CP.Yr2.mod4.se<- se(N.Sub2.Spp2$delta.15.N)

CP.Yr3.mod5<- mean(C.Sub3.Spp2$delta.13.C)
CP.Yr3.mod5.se<- se(C.Sub3.Spp2$delta.13.C)

CP.Yr3.mod6<- mean(N.Sub3.Spp2$delta.15.N)
CP.Yr3.mod6.se<- se(N.Sub3.Spp2$delta.15.N)

#-----
# GEPE

GP.Yr1.mod1<- mean(C.Sub1.Spp3$delta.13.C)
GP.Yr1.mod1.se<- se(C.Sub1.Spp3$delta.13.C)

GP.Yr1.mod2<- mean(N.Sub1.Spp3$delta.15.N)
GP.Yr1.mod2.se<- se(N.Sub1.Spp3$delta.15.N)

GP.Yr2.mod3<- mean(C.Sub2.Spp3$delta.13.C)
GP.Yr2.mod3.se<- se(C.Sub2.Spp3$delta.13.C)

GP.Yr2.mod4<- mean(N.Sub2.Spp3$delta.15.N)
GP.Yr2.mod4.se<- se(N.Sub2.Spp3$delta.15.N)

GP.Yr3.mod5<- mean(C.Sub3.Spp3$delta.13.C)
GP.Yr3.mod5.se<- se(C.Sub3.Spp3$delta.13.C)

GP.Yr3.mod6<- mean(N.Sub3.Spp3$delta.15.N)
GP.Yr3.mod6.se<- se(N.Sub3.Spp3$delta.15.N)

#-----
# ADPE. Plot each year by isotopes

AP.x1<- AP.Yr1.mod1
AP.y1<- AP.Yr1.mod2
AP.x1.se<- AP.Yr1.mod1.se
AP.y1.se<- AP.Yr1.mod2.se

AP.x2<- AP.Yr2.mod3
AP.y2<- AP.Yr2.mod4
AP.x2.se<- AP.Yr2.mod3.se
AP.y2.se<- AP.Yr2.mod4.se

AP.x3<- AP.Yr3.mod5
AP.y3<- AP.Yr3.mod6
AP.x3.se<- AP.Yr3.mod5.se
AP.y3.se<- AP.Yr3.mod6.se

#-----
# CHPE. Plot each year by isotopes

CP.x1<- CP.Yr1.mod1
CP.y1<- CP.Yr1.mod2
CP.x1.se<- CP.Yr1.mod1.se
CP.y1.se<- CP.Yr1.mod2.se

CP.x2<- CP.Yr2.mod3
CP.y2<- CP.Yr2.mod4
CP.x2.se<- CP.Yr2.mod3.se
CP.y2.se<- CP.Yr2.mod4.se

CP.x3<- CP.Yr3.mod5
CP.y3<- CP.Yr3.mod6
CP.x3.se<- CP.Yr3.mod5.se
CP.y3.se<- CP.Yr3.mod6.se

#-----
# GEPE. Plot each year by isotopes

GP.x1<- GP.Yr1.mod1
GP.y1<- GP.Yr1.mod2
GP.x1.se<- GP.Yr1.mod1.se
GP.y1.se<- GP.Yr1.mod2.se

GP.x2<- GP.Yr2.mod3
GP.y2<- GP.Yr2.mod4
GP.x2.se<- GP.Yr2.mod3.se
GP.y2.se<- GP.Yr2.mod4.se

GP.x3<- GP.Yr3.mod5
GP.y3<- GP.Yr3.mod6
GP.x3.se<- GP.Yr3.mod5.se
GP.y3.se<- GP.Yr3.mod6.se

#-----
# Plot

tiff(file="CKs.D5.CN.PAL.Year2.tiff",res=300,width=15,height=15,unit="cm")
par(mfrow=c(1,3))

plot1<- plot(AP.y1~AP.x1, main="2007/08", cex.main=1.5, pch=19, col="black", xlim=c(-28.0,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(AP.y1~AP.x1)
arrows(AP.x1, AP.y1-AP.y1.se, AP.x1, AP.y1+AP.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x1-AP.x1.se, AP.y1, AP.x1+AP.x1.se, AP.y1, code=3, angle=90, length=0.1, col="black")

points(CP.y1~CP.x1, pch=19, col="orangered")
arrows(CP.x1, CP.y1-CP.y1.se, CP.x1, CP.y1+CP.y1.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x1-CP.x1.se, CP.y1, CP.x1+CP.x1.se, CP.y1, code=3, angle=90, length=0.1, col="orangered")

points(GP.y1~GP.x1, pch=19, col="darkblue")
arrows(GP.x1, GP.y1-GP.y1.se, GP.x1, GP.y1+GP.y1.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x1-GP.x1.se, GP.y1, GP.x1+GP.x1.se, GP.y1, code=3, angle=90, length=0.1, col="darkblue")

plot2<- plot(AP.y2~AP.x2, main="2008/09", cex.main=1.5, pch=19, col="black", xlim=c(-28.0,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(AP.y2~AP.x2)
arrows(AP.x2, AP.y2-AP.y2.se, AP.x2, AP.y2+AP.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x2-AP.x2.se, AP.y2, AP.x2+AP.x2.se, AP.y2, code=3, angle=90, length=0.1, col="black")

points(CP.y2~CP.x2, pch=19, col="orangered")
arrows(CP.x2, CP.y2-CP.y2.se, CP.x2, CP.y2+CP.y2.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x2-CP.x2.se, CP.y2, CP.x2+CP.x2.se, CP.y2, code=3, angle=90, length=0.1, col="orangered")

points(GP.y2~GP.x2, pch=19, col="darkblue")
arrows(GP.x2, GP.y2-GP.y2.se, GP.x2, GP.y2+GP.y2.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x2-GP.x2.se, GP.y2, GP.x2+GP.x2.se, GP.y2, code=3, angle=90, length=0.1, col="darkblue")

plot3<- plot(AP.y3~AP.x3, main="2009/10", cex.main=1.5, pch=19, col="black", xlim=c(-28.0,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(AP.y3~AP.x3)
arrows(AP.x3, AP.y3-AP.y3.se, AP.x3, AP.y3+AP.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x3-AP.x3.se, AP.y3, AP.x3+AP.x3.se, AP.y3, code=3, angle=90, length=0.1, col="black")

points(CP.y3~CP.x3, pch=19, col="orangered")
arrows(CP.x3, CP.y3-CP.y3.se, CP.x3, CP.y3+CP.y3.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x3-CP.x3.se, CP.y3, CP.x3+CP.x3.se, CP.y3, code=3, angle=90, length=0.1, col="orangered")

points(GP.y3~GP.x3, pch=19, col="darkblue")
arrows(GP.x3, GP.y3-GP.y3.se, GP.x3, GP.y3+GP.y3.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x3-GP.x3.se, GP.y3, GP.x3+GP.x3.se, GP.y3, code=3, angle=90, length=0.1, col="darkblue")

dev.off()

#-----
# Species isotope graphs by year, Cks.D15.PAL

#-----
# Reset R, clear all objects
rm(list=ls()) 

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

DataSet<- read.csv("Cks.D15.PAL.INS.DataSet.Sub2_FINAL.csv")

head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

DataSet$Year.2<- as.factor(DataSet$Year.2)
DataSet$Nest<- as.factor(DataSet$Nest)
str(DataSet)
DataSet

range(DataSet$delta.15.N)
range(DataSet$delta.13.C)

#-----
DataSet1<- DataSet

head(DataSet1)
summary(DataSet1)
nrow(DataSet1)
str(DataSet1)

#-----
# Subset DataSets into spp and year, by isotope. Can use same dataset for both N and C.

#-----
# DataSet1.N

N.Sub1<- subset(DataSet1,Year.2=="1")
N.Sub2<- subset(DataSet1,Year.2=="2")
N.Sub3<- subset(DataSet1,Year.2=="3") 

N.Sub1.Spp1<- subset(N.Sub1,Spp=="ADPE")
N.Sub1.Spp2<- subset(N.Sub1,Spp=="CHPE")
N.Sub1.Spp3<- subset(N.Sub1,Spp=="GEPE")

N.Sub2.Spp1<- subset(N.Sub2,Spp=="ADPE")
N.Sub2.Spp2<- subset(N.Sub2,Spp=="CHPE")
N.Sub2.Spp3<- subset(N.Sub2,Spp=="GEPE")

N.Sub3.Spp1<- subset(N.Sub3,Spp=="ADPE")
N.Sub3.Spp2<- subset(N.Sub3,Spp=="CHPE")
N.Sub3.Spp3<- subset(N.Sub3,Spp=="GEPE")

#-----
# DataSet1.C

C.Sub1<- subset(DataSet1,Year.2=="1")
C.Sub2<- subset(DataSet1,Year.2=="2")
C.Sub3<- subset(DataSet1,Year.2=="3") 

C.Sub1.Spp1<- subset(C.Sub1,Spp=="ADPE")
C.Sub1.Spp2<- subset(C.Sub1,Spp=="CHPE")
C.Sub1.Spp3<- subset(C.Sub1,Spp=="GEPE")

C.Sub2.Spp1<- subset(C.Sub2,Spp=="ADPE")
C.Sub2.Spp2<- subset(C.Sub2,Spp=="CHPE")
C.Sub2.Spp3<- subset(C.Sub2,Spp=="GEPE")

C.Sub3.Spp1<- subset(C.Sub3,Spp=="ADPE")
C.Sub3.Spp2<- subset(C.Sub3,Spp=="CHPE")
C.Sub3.Spp3<- subset(C.Sub3,Spp=="GEPE")

#-----
# Averages for each spp by year for each isotope, and standard error. The variance estimate (var) is a built in function in R. The stdev is the sqrt of the variance. The se function was checked in excel and gave the same value.

variance<- function(x) sum((x-mean(x))^2)/(length(x)-1)
stdev<- function(x) sqrt(sum((x-mean(x))^2)/(length(x)-1))
se<- function(x) sqrt(var(x)/length(x))

# ADPE

AP.Yr1.mod1<- mean(C.Sub1.Spp1$delta.13.C)
AP.Yr1.mod1.se<- se(C.Sub1.Spp1$delta.13.C)

AP.Yr1.mod2<- mean(N.Sub1.Spp1$delta.15.N)
AP.Yr1.mod2.se<- se(N.Sub1.Spp1$delta.15.N)

AP.Yr2.mod3<- mean(C.Sub2.Spp1$delta.13.C)
AP.Yr2.mod3.se<- se(C.Sub2.Spp1$delta.13.C)

AP.Yr2.mod4<- mean(N.Sub2.Spp1$delta.15.N)
AP.Yr2.mod4.se<- se(N.Sub2.Spp1$delta.15.N)

AP.Yr3.mod5<- mean(C.Sub3.Spp1$delta.13.C)
AP.Yr3.mod5.se<- se(C.Sub3.Spp1$delta.13.C)

AP.Yr3.mod6<- mean(N.Sub3.Spp1$delta.15.N)
AP.Yr3.mod6.se<- se(N.Sub3.Spp1$delta.15.N)

#-----
# CHPE

CP.Yr1.mod1<- mean(C.Sub1.Spp2$delta.13.C)
CP.Yr1.mod1.se<- se(C.Sub1.Spp2$delta.13.C)

CP.Yr1.mod2<- mean(N.Sub1.Spp2$delta.15.N)
CP.Yr1.mod2.se<- se(N.Sub1.Spp2$delta.15.N)

CP.Yr2.mod3<- mean(C.Sub2.Spp2$delta.13.C)
CP.Yr2.mod3.se<- se(C.Sub2.Spp2$delta.13.C)

CP.Yr2.mod4<- mean(N.Sub2.Spp2$delta.15.N)
CP.Yr2.mod4.se<- se(N.Sub2.Spp2$delta.15.N)

CP.Yr3.mod5<- mean(C.Sub3.Spp2$delta.13.C)
CP.Yr3.mod5.se<- se(C.Sub3.Spp2$delta.13.C)

CP.Yr3.mod6<- mean(N.Sub3.Spp2$delta.15.N)
CP.Yr3.mod6.se<- se(N.Sub3.Spp2$delta.15.N)

#-----
# GEPE

GP.Yr1.mod1<- mean(C.Sub1.Spp3$delta.13.C)
GP.Yr1.mod1.se<- se(C.Sub1.Spp3$delta.13.C)

GP.Yr1.mod2<- mean(N.Sub1.Spp3$delta.15.N)
GP.Yr1.mod2.se<- se(N.Sub1.Spp3$delta.15.N)

GP.Yr2.mod3<- mean(C.Sub2.Spp3$delta.13.C)
GP.Yr2.mod3.se<- se(C.Sub2.Spp3$delta.13.C)

GP.Yr2.mod4<- mean(N.Sub2.Spp3$delta.15.N)
GP.Yr2.mod4.se<- se(N.Sub2.Spp3$delta.15.N)

GP.Yr3.mod5<- mean(C.Sub3.Spp3$delta.13.C)
GP.Yr3.mod5.se<- se(C.Sub3.Spp3$delta.13.C)

GP.Yr3.mod6<- mean(N.Sub3.Spp3$delta.15.N)
GP.Yr3.mod6.se<- se(N.Sub3.Spp3$delta.15.N)

#-----
# ADPE. Plot each year by isotopes

AP.x1<- AP.Yr1.mod1
AP.y1<- AP.Yr1.mod2
AP.x1.se<- AP.Yr1.mod1.se
AP.y1.se<- AP.Yr1.mod2.se

AP.x2<- AP.Yr2.mod3
AP.y2<- AP.Yr2.mod4
AP.x2.se<- AP.Yr2.mod3.se
AP.y2.se<- AP.Yr2.mod4.se

AP.x3<- AP.Yr3.mod5
AP.y3<- AP.Yr3.mod6
AP.x3.se<- AP.Yr3.mod5.se
AP.y3.se<- AP.Yr3.mod6.se

#-----
# CHPE. Plot each year by isotopes

CP.x1<- CP.Yr1.mod1
CP.y1<- CP.Yr1.mod2
CP.x1.se<- CP.Yr1.mod1.se
CP.y1.se<- CP.Yr1.mod2.se

CP.x2<- CP.Yr2.mod3
CP.y2<- CP.Yr2.mod4
CP.x2.se<- CP.Yr2.mod3.se
CP.y2.se<- CP.Yr2.mod4.se

CP.x3<- CP.Yr3.mod5
CP.y3<- CP.Yr3.mod6
CP.x3.se<- CP.Yr3.mod5.se
CP.y3.se<- CP.Yr3.mod6.se

#-----
# GEPE. Plot each year by isotopes

GP.x1<- GP.Yr1.mod1
GP.y1<- GP.Yr1.mod2
GP.x1.se<- GP.Yr1.mod1.se
GP.y1.se<- GP.Yr1.mod2.se

GP.x2<- GP.Yr2.mod3
GP.y2<- GP.Yr2.mod4
GP.x2.se<- GP.Yr2.mod3.se
GP.y2.se<- GP.Yr2.mod4.se

GP.x3<- GP.Yr3.mod5
GP.y3<- GP.Yr3.mod6
GP.x3.se<- GP.Yr3.mod5.se
GP.y3.se<- GP.Yr3.mod6.se


#-----
# Plot

tiff(file="CKs.D15.CN.PAL.Year2.tiff",res=300,width=15,height=15,unit="cm")
par(mfrow=c(1,3))

plot1<- plot(AP.y1~AP.x1, main="2007/08", cex.main=1.5, pch=19, col="black", xlim=c(-28.0,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.5,5.5,4.1,2.1)))
points(AP.y1~AP.x1)
arrows(AP.x1, AP.y1-AP.y1.se, AP.x1, AP.y1+AP.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x1-AP.x1.se, AP.y1, AP.x1+AP.x1.se, AP.y1, code=3, angle=90, length=0.1, col="black")

points(CP.y1~CP.x1, pch=19, col="orangered")
arrows(CP.x1, CP.y1-CP.y1.se, CP.x1, CP.y1+CP.y1.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x1-CP.x1.se, CP.y1, CP.x1+CP.x1.se, CP.y1, code=3, angle=90, length=0.1, col="orangered")

points(GP.y1~GP.x1, pch=19, col="darkblue")
arrows(GP.x1, GP.y1-GP.y1.se, GP.x1, GP.y1+GP.y1.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x1-GP.x1.se, GP.y1, GP.x1+GP.x1.se, GP.y1, code=3, angle=90, length=0.1, col="darkblue")

plot2<- plot(AP.y2~AP.x2, main="2008/09", cex.main=1.5, pch=19, col="black", xlim=c(-28.0,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.5,5.5,4.1,2.1)))
points(AP.y2~AP.x2)
arrows(AP.x2, AP.y2-AP.y2.se, AP.x2, AP.y2+AP.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x2-AP.x2.se, AP.y2, AP.x2+AP.x2.se, AP.y2, code=3, angle=90, length=0.1, col="black")

points(CP.y2~CP.x2, pch=19, col="orangered")
arrows(CP.x2, CP.y2-CP.y2.se, CP.x2, CP.y2+CP.y2.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x2-CP.x2.se, CP.y2, CP.x2+CP.x2.se, CP.y2, code=3, angle=90, length=0.1, col="orangered")

points(GP.y2~GP.x2, pch=19, col="darkblue")
arrows(GP.x2, GP.y2-GP.y2.se, GP.x2, GP.y2+GP.y2.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x2-GP.x2.se, GP.y2, GP.x2+GP.x2.se, GP.y2, code=3, angle=90, length=0.1, col="darkblue")

plot3<- plot(AP.y3~AP.x3, main="2009/10", cex.main=1.5, pch=19, col="black", xlim=c(-28.0,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.5,5.5,4.1,2.1)))
points(AP.y3~AP.x3)
arrows(AP.x3, AP.y3-AP.y3.se, AP.x3, AP.y3+AP.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x3-AP.x3.se, AP.y3, AP.x3+AP.x3.se, AP.y3, code=3, angle=90, length=0.1, col="black")

points(CP.y3~CP.x3, pch=19, col="orangered")
arrows(CP.x3, CP.y3-CP.y3.se, CP.x3, CP.y3+CP.y3.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x3-CP.x3.se, CP.y3, CP.x3+CP.x3.se, CP.y3, code=3, angle=90, length=0.1, col="orangered")

points(GP.y3~GP.x3, pch=19, col="darkblue")
arrows(GP.x3, GP.y3-GP.y3.se, GP.x3, GP.y3+GP.y3.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x3-GP.x3.se, GP.y3, GP.x3+GP.x3.se, GP.y3, code=3, angle=90, length=0.1, col="darkblue")

dev.off()

#-----
# Species isotope graphs by year. Cks.5Wk.PAL

#-----
# Reset R, clear all objects
rm(list=ls()) 

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

DataSet<- read.csv("Cks.5wk.PAL.INS.DataSet.Sub7_FINAL.csv")

head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

range(DataSet$delta.15.N)
range(DataSet$delta.13.C)

#-----
DataSet1<- DataSet

head(DataSet1)
summary(DataSet1)
nrow(DataSet1)
str(DataSet1)

#-----
# Subset DataSets into spp and year, by isotope

#-----
# DataSet1.N

N.Sub1<- subset(DataSet1,Year.2=="1")
N.Sub2<- subset(DataSet1,Year.2=="2")
N.Sub3<- subset(DataSet1,Year.2=="3") 

N.Sub1.Spp1<- subset(N.Sub1,Spp=="ADPE")
N.Sub1.Spp2<- subset(N.Sub1,Spp=="CHPE")
N.Sub1.Spp3<- subset(N.Sub1,Spp=="GEPE")

N.Sub2.Spp1<- subset(N.Sub2,Spp=="ADPE")
N.Sub2.Spp2<- subset(N.Sub2,Spp=="CHPE")
N.Sub2.Spp3<- subset(N.Sub2,Spp=="GEPE")

N.Sub3.Spp1<- subset(N.Sub3,Spp=="ADPE")
N.Sub3.Spp2<- subset(N.Sub3,Spp=="CHPE")
N.Sub3.Spp3<- subset(N.Sub3,Spp=="GEPE")

#-----
# DataSet1.C

C.Sub1<- subset(DataSet1,Year.2=="1")
C.Sub2<- subset(DataSet1,Year.2=="2")
C.Sub3<- subset(DataSet1,Year.2=="3") 

C.Sub1.Spp1<- subset(C.Sub1,Spp=="ADPE")
C.Sub1.Spp2<- subset(C.Sub1,Spp=="CHPE")
C.Sub1.Spp3<- subset(C.Sub1,Spp=="GEPE")

C.Sub2.Spp1<- subset(C.Sub2,Spp=="ADPE")
C.Sub2.Spp2<- subset(C.Sub2,Spp=="CHPE")
C.Sub2.Spp3<- subset(C.Sub2,Spp=="GEPE")

C.Sub3.Spp1<- subset(C.Sub3,Spp=="ADPE")
C.Sub3.Spp2<- subset(C.Sub3,Spp=="CHPE")
C.Sub3.Spp3<- subset(C.Sub3,Spp=="GEPE")

#-----
# Averages for each spp by year for each isotope, and standard error. The variance estimate (var) is a built in function in R. The stdev is the sqrt of the variance. The se function was checked in excel and gave the same value.

variance<- function(x) sum((x-mean(x))^2)/(length(x)-1)
stdev<- function(x) sqrt(sum((x-mean(x))^2)/(length(x)-1))
se<- function(x) sqrt(var(x)/length(x))

# ADPE

AP.Yr1.mod1<- mean(C.Sub1.Spp1$delta.13.C)
AP.Yr1.mod1.se<- se(C.Sub1.Spp1$delta.13.C)

AP.Yr1.mod2<- mean(N.Sub1.Spp1$delta.15.N)
AP.Yr1.mod2.se<- se(N.Sub1.Spp1$delta.15.N)

AP.Yr2.mod3<- mean(C.Sub2.Spp1$delta.13.C)
AP.Yr2.mod3.se<- se(C.Sub2.Spp1$delta.13.C)

AP.Yr2.mod4<- mean(N.Sub2.Spp1$delta.15.N)
AP.Yr2.mod4.se<- se(N.Sub2.Spp1$delta.15.N)

AP.Yr3.mod5<- mean(C.Sub3.Spp1$delta.13.C)
AP.Yr3.mod5.se<- se(C.Sub3.Spp1$delta.13.C)

AP.Yr3.mod6<- mean(N.Sub3.Spp1$delta.15.N)
AP.Yr3.mod6.se<- se(N.Sub3.Spp1$delta.15.N)

#-----
# CHPE

CP.Yr1.mod1<- mean(C.Sub1.Spp2$delta.13.C)
CP.Yr1.mod1.se<- se(C.Sub1.Spp2$delta.13.C)

CP.Yr1.mod2<- mean(N.Sub1.Spp2$delta.15.N)
CP.Yr1.mod2.se<- se(N.Sub1.Spp2$delta.15.N)

CP.Yr2.mod3<- mean(C.Sub2.Spp2$delta.13.C)
CP.Yr2.mod3.se<- se(C.Sub2.Spp2$delta.13.C)

CP.Yr2.mod4<- mean(N.Sub2.Spp2$delta.15.N)
CP.Yr2.mod4.se<- se(N.Sub2.Spp2$delta.15.N)

CP.Yr3.mod5<- mean(C.Sub3.Spp2$delta.13.C)
CP.Yr3.mod5.se<- se(C.Sub3.Spp2$delta.13.C)

CP.Yr3.mod6<- mean(N.Sub3.Spp2$delta.15.N)
CP.Yr3.mod6.se<- se(N.Sub3.Spp2$delta.15.N)

#-----
# GEPE

GP.Yr1.mod1<- mean(C.Sub1.Spp3$delta.13.C)
GP.Yr1.mod1.se<- se(C.Sub1.Spp3$delta.13.C)

GP.Yr1.mod2<- mean(N.Sub1.Spp3$delta.15.N)
GP.Yr1.mod2.se<- se(N.Sub1.Spp3$delta.15.N)

GP.Yr2.mod3<- mean(C.Sub2.Spp3$delta.13.C)
GP.Yr2.mod3.se<- se(C.Sub2.Spp3$delta.13.C)

GP.Yr2.mod4<- mean(N.Sub2.Spp3$delta.15.N)
GP.Yr2.mod4.se<- se(N.Sub2.Spp3$delta.15.N)

GP.Yr3.mod5<- mean(C.Sub3.Spp3$delta.13.C)
GP.Yr3.mod5.se<- se(C.Sub3.Spp3$delta.13.C)

GP.Yr3.mod6<- mean(N.Sub3.Spp3$delta.15.N)
GP.Yr3.mod6.se<- se(N.Sub3.Spp3$delta.15.N)

#-----
# ADPE. Plot each year by isotopes

AP.x1<- AP.Yr1.mod1
AP.y1<- AP.Yr1.mod2
AP.x1.se<- AP.Yr1.mod1.se
AP.y1.se<- AP.Yr1.mod2.se

AP.x2<- AP.Yr2.mod3
AP.y2<- AP.Yr2.mod4
AP.x2.se<- AP.Yr2.mod3.se
AP.y2.se<- AP.Yr2.mod4.se

AP.x3<- AP.Yr3.mod5
AP.y3<- AP.Yr3.mod6
AP.x3.se<- AP.Yr3.mod5.se
AP.y3.se<- AP.Yr3.mod6.se

#-----
# CHPE. Plot each year by isotopes

CP.x1<- CP.Yr1.mod1
CP.y1<- CP.Yr1.mod2
CP.x1.se<- CP.Yr1.mod1.se
CP.y1.se<- CP.Yr1.mod2.se

CP.x2<- CP.Yr2.mod3
CP.y2<- CP.Yr2.mod4
CP.x2.se<- CP.Yr2.mod3.se
CP.y2.se<- CP.Yr2.mod4.se

CP.x3<- CP.Yr3.mod5
CP.y3<- CP.Yr3.mod6
CP.x3.se<- CP.Yr3.mod5.se
CP.y3.se<- CP.Yr3.mod6.se

#-----
# GEPE. Plot each year by isotopes

GP.x1<- GP.Yr1.mod1
GP.y1<- GP.Yr1.mod2
GP.x1.se<- GP.Yr1.mod1.se
GP.y1.se<- GP.Yr1.mod2.se

GP.x2<- GP.Yr2.mod3
GP.y2<- GP.Yr2.mod4
GP.x2.se<- GP.Yr2.mod3.se
GP.y2.se<- GP.Yr2.mod4.se

GP.x3<- GP.Yr3.mod5
GP.y3<- GP.Yr3.mod6
GP.x3.se<- GP.Yr3.mod5.se
GP.y3.se<- GP.Yr3.mod6.se

#-----
# Plot

tiff(file="CKs.5wk.CN.PAL.Year2.tiff",res=300,width=15,height=15,unit="cm")
par(mfrow=c(1,3))

plot1<- plot(AP.y1~AP.x1, main="2007/08", cex.main=1.5, pch=19, col="black", xlim=c(-28,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y1~CP.x1, pch=19, col="orangered1")
points(GP.y1~GP.x1, pch=19, col="darkblue")

points(AP.x1, AP.y1, pch=23, col="black", bg="black")
arrows(AP.x1, AP.y1-AP.y1.se, AP.x1, AP.y1+AP.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x1-AP.x1.se, AP.y1, AP.x1+AP.x1.se, AP.y1, code=3, angle=90, length=0.1, col="black")

points(CP.x1, CP.y1, pch=23, col="orangered", bg="orangered")
arrows(CP.x1, CP.y1-CP.y1.se, CP.x1, CP.y1+CP.y1.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x1-CP.x1.se, CP.y1, CP.x1+CP.x1.se, CP.y1, code=3, angle=90, length=0.1, col="orangered")

points(GP.x1, GP.y1, pch=23, col="darkblue", bg="darkblue")
arrows(GP.x1, GP.y1-GP.y1.se, GP.x1, GP.y1+GP.y1.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x1-GP.x1.se, GP.y1, GP.x1+GP.x1.se, GP.y1, code=3, angle=90, length=0.1, col="darkblue")

plot2<- plot(AP.y2~AP.x2, main="2008/09", cex.main=1.5, pch=19, col="black", xlim=c(-28,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y2~CP.x2, pch=19, col="orangered1")
points(GP.y2~GP.x2, pch=19, col="darkblue")

points(AP.x2, AP.y2, pch=23, col="black", bg="black")
arrows(AP.x2, AP.y2-AP.y2.se, AP.x2, AP.y2+AP.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x2-AP.x2.se, AP.y2, AP.x2+AP.x2.se, AP.y2, code=3, angle=90, length=0.1, col="black")

points(CP.x2, CP.y2, pch=23, col="orangered", bg="orangered")
arrows(CP.x2, CP.y2-CP.y2.se, CP.x2, CP.y2+CP.y2.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x2-CP.x2.se, CP.y2, CP.x2+CP.x2.se, CP.y2, code=3, angle=90, length=0.1, col="orangered")

points(GP.x2, GP.y2, pch=23, col="darkblue", bg="darkblue")
arrows(GP.x2, GP.y2-GP.y2.se, GP.x2, GP.y2+GP.y2.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x2-GP.x2.se, GP.y2, GP.x2+GP.x2.se, GP.y2, code=3, angle=90, length=0.1, col="darkblue")

plot3<- plot(AP.y3~AP.x3, main="2009/10", cex.main=1.5, pch=19, col="black", xlim=c(-28,-24), ylim=c(6.8,9.6), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y3~CP.x3, pch=19, col="orangered1")
points(GP.y3~GP.x3, pch=19, col="darkblue")

points(AP.x3, AP.y3, pch=23, col="black", bg="black")
arrows(AP.x3, AP.y3-AP.y3.se, AP.x3, AP.y3+AP.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x3-AP.x3.se, AP.y3, AP.x3+AP.x3.se, AP.y3, code=3, angle=90, length=0.1, col="black")

points(CP.x3, CP.y3, pch=23, col="orangered", bg="orangered")
arrows(CP.x3, CP.y3-CP.y3.se, CP.x3, CP.y3+CP.y3.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x3-CP.x3.se, CP.y3, CP.x3+CP.x3.se, CP.y3, code=3, angle=90, length=0.1, col="orangered")

points(GP.x3, GP.y3, pch=23, col="darkblue", bg="darkblue")
arrows(GP.x3, GP.y3-GP.y3.se, GP.x3, GP.y3+GP.y3.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x3-GP.x3.se, GP.y3, GP.x3+GP.x3.se, GP.y3, code=3, angle=90, length=0.1, col="darkblue")

dev.off()

#-----
# Species isotope graphs by year. Cks.5Wk.WAP

#-----
# Reset R, clear all objects
rm(list=ls()) 

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

DataSet<- read.csv("Cks.5wk.WAP.INS.DataSet.Sub3_FINAL.csv")

head(DataSet)
summary(DataSet)
nrow(DataSet)
str(DataSet)

DataSet$Year.2<- as.factor(DataSet$Year.2)
str(DataSet)
DataSet

range(DataSet$delta.15.N)

#-----
DataSet1<- DataSet

head(DataSet1)
summary(DataSet1)
nrow(DataSet1)
str(DataSet1)

#-----
# Subset DataSets into spp and year, by isotope

#-----
# DataSet1.N

N.Sub1<- subset(DataSet1,Year.2=="1")
N.Sub2<- subset(DataSet1,Year.2=="2")
N.Sub3<- subset(DataSet1,Year.2=="3") 

N.Sub1.Loc1<- subset(N.Sub1,Location=="PAL")
N.Sub1.Loc2<- subset(N.Sub1,Location=="AVI")

N.Sub2.Loc1<- subset(N.Sub2,Location=="PAL")
N.Sub2.Loc2<- subset(N.Sub2,Location=="AVI")

N.Sub3.Loc1<- subset(N.Sub3,Location=="PAL")
N.Sub3.Loc2<- subset(N.Sub3,Location=="AVI")
N.Sub3.Loc3<- subset(N.Sub3,Location=="CHA")

#-----
# DataSet1.C

list.files()
DataSet2<- read.csv("Cks.5wk.WAP.INS.DataSet.Sub4_FINAL.csv")

head(DataSet2)
summary(DataSet2)
nrow(DataSet2)
str(DataSet2)

DataSet2$Year.2<- as.factor(DataSet2$Year.2)
str(DataSet2)
DataSet2

range(DataSet$delta.13.C)

C.Sub1<- subset(DataSet2,Year.2=="1")
C.Sub2<- subset(DataSet2,Year.2=="2")
C.Sub3<- subset(DataSet2,Year.2=="3") 

C.Sub1.Loc1<- subset(C.Sub1,Location=="PAL")
C.Sub1.Loc2<- subset(C.Sub1,Location=="AVI")

C.Sub2.Loc1<- subset(C.Sub2,Location=="PAL")
C.Sub2.Loc2<- subset(C.Sub2,Location=="AVI")

C.Sub3.Loc1<- subset(C.Sub3,Location=="PAL")
C.Sub3.Loc2<- subset(C.Sub3,Location=="AVI")
C.Sub3.Loc3<- subset(C.Sub3,Location=="CHA")

#-----
# Averages for each spp by year for each isotope, and standard error. The variance estimate (var) is a built in function in R. The stdev is the sqrt of the variance. The se function was checked in excel and gave the same value.

variance<- function(x) sum((x-mean(x))^2)/(length(x)-1)
stdev<- function(x) sqrt(sum((x-mean(x))^2)/(length(x)-1))
se<- function(x) sqrt(var(x)/length(x))

# PAL

PAL.Yr1.mod1<- mean(C.Sub1.Loc1$delta.13.C)
PAL.Yr1.mod1.se<- se(C.Sub1.Loc1$delta.13.C)

PAL.Yr1.mod2<- mean(N.Sub1.Loc1$delta.15.N)
PAL.Yr1.mod2.se<- se(N.Sub1.Loc1$delta.15.N)

PAL.Yr2.mod3<- mean(C.Sub2.Loc1$delta.13.C)
PAL.Yr2.mod3.se<- se(C.Sub2.Loc1$delta.13.C)

PAL.Yr2.mod4<- mean(N.Sub2.Loc1$delta.15.N)
PAL.Yr2.mod4.se<- se(N.Sub2.Loc1$delta.15.N)

PAL.Yr3.mod5<- mean(C.Sub3.Loc1$delta.13.C)
PAL.Yr3.mod5.se<- se(C.Sub3.Loc1$delta.13.C)

PAL.Yr3.mod6<- mean(N.Sub3.Loc1$delta.15.N)
PAL.Yr3.mod6.se<- se(N.Sub3.Loc1$delta.15.N)

#-----
# AVI

AVI.Yr1.mod1<- mean(C.Sub1.Loc2$delta.13.C)
AVI.Yr1.mod1.se<- se(C.Sub1.Loc2$delta.13.C)

AVI.Yr1.mod2<- mean(N.Sub1.Loc2$delta.15.N)
AVI.Yr1.mod2.se<- se(N.Sub1.Loc2$delta.15.N)

AVI.Yr2.mod3<- mean(C.Sub2.Loc2$delta.13.C)
AVI.Yr2.mod3.se<- se(C.Sub2.Loc2$delta.13.C)

AVI.Yr2.mod4<- mean(N.Sub2.Loc2$delta.15.N)
AVI.Yr2.mod4.se<- se(N.Sub2.Loc2$delta.15.N)

AVI.Yr3.mod5<- mean(C.Sub3.Loc2$delta.13.C)
AVI.Yr3.mod5.se<- se(C.Sub3.Loc2$delta.13.C)

AVI.Yr3.mod6<- mean(N.Sub3.Loc2$delta.15.N)
AVI.Yr3.mod6.se<- se(N.Sub3.Loc2$delta.15.N)

#-----
# CHA

CHA.Yr3.mod5<- mean(C.Sub3.Loc3$delta.13.C)
CHA.Yr3.mod5.se<- se(C.Sub3.Loc3$delta.13.C)

CHA.Yr3.mod6<- mean(N.Sub3.Loc3$delta.15.N)
CHA.Yr3.mod6.se<- se(N.Sub3.Loc3$delta.15.N)

#-----
# PAL. Plot each location by isotopes

PAL.x1<- PAL.Yr1.mod1
PAL.y1<- PAL.Yr1.mod2
PAL.x1.se<- PAL.Yr1.mod1.se
PAL.y1.se<- PAL.Yr1.mod2.se

PAL.x2<- PAL.Yr2.mod3
PAL.y2<- PAL.Yr2.mod4
PAL.x2.se<- PAL.Yr2.mod3.se
PAL.y2.se<- PAL.Yr2.mod4.se

PAL.x3<- PAL.Yr3.mod5
PAL.y3<- PAL.Yr3.mod6
PAL.x3.se<- PAL.Yr3.mod5.se
PAL.y3.se<- PAL.Yr3.mod6.se

#-----
# AVI. Plot each location by isotopes

AVI.x1<- AVI.Yr1.mod1
AVI.y1<- AVI.Yr1.mod2
AVI.x1.se<- AVI.Yr1.mod1.se
AVI.y1.se<- AVI.Yr1.mod2.se

AVI.x2<- AVI.Yr2.mod3
AVI.y2<- AVI.Yr2.mod4
AVI.x2.se<- AVI.Yr2.mod3.se
AVI.y2.se<- AVI.Yr2.mod4.se

AVI.x3<- AVI.Yr3.mod5
AVI.y3<- AVI.Yr3.mod6
AVI.x3.se<- AVI.Yr3.mod5.se
AVI.y3.se<- AVI.Yr3.mod6.se

#-----
# CHA. Plot each location by isotopes

CHA.x3<- CHA.Yr3.mod5
CHA.y3<- CHA.Yr3.mod6
CHA.x3.se<- CHA.Yr3.mod5.se
CHA.y3.se<- CHA.Yr3.mod6.se

#-----
# Plot

tiff(file="CKs.CN.Year.Loc.tiff",res=150,width=15,height=15,unit="cm")
par(mfrow=c(1,3))

plot1<- plot(PAL.y1~PAL.x1, main="2007/08", cex.main=1.5, pch=21, cex=1.4, col="black", bg="black", xlim=c(-28.5,-22.5), ylim=c(7.0,9.2), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(AVI.y1~AVI.x1, pch=22, cex=1.4, col="grey", bg="grey")

arrows(PAL.x1, PAL.y1-PAL.y1.se, PAL.x1, PAL.y1+PAL.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(PAL.x1-PAL.x1.se, PAL.y1, PAL.x1+PAL.x1.se, PAL.y1, code=3, angle=90, length=0.1, col="black")
points(PAL.x1, PAL.y1, pch=21, cex=1.4, col="black", bg="black")

arrows(AVI.x1, AVI.y1-AVI.y1.se, AVI.x1, AVI.y1+AVI.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(AVI.x1-AVI.x1.se, AVI.y1, AVI.x1+AVI.x1.se, AVI.y1, code=3, angle=90, length=0.1, col="black")
points(AVI.x1, AVI.y1, pch=22, cex=1.4, col="black", bg="grey")


plot2<- plot(PAL.y2~PAL.x2, main="2008/09", cex.main=1.5, pch=21, cex=1.4, col="black", bg="black", xlim=c(-28.5,-22.5), ylim=c(7.0,9.2), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(AVI.y2~AVI.x2, pch=22, cex=1.4, col="grey", bg="grey")

arrows(PAL.x2, PAL.y2-PAL.y2.se, PAL.x2, PAL.y2+PAL.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(PAL.x2-PAL.x2.se, PAL.y2, PAL.x2+PAL.x2.se, PAL.y2, code=3, angle=90, length=0.1, col="black")
points(PAL.x2, PAL.y2, pch=21, cex=1.4, col="black", bg="black")

arrows(AVI.x2, AVI.y2-AVI.y2.se, AVI.x2, AVI.y2+AVI.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(AVI.x2-AVI.x2.se, AVI.y2, AVI.x2+AVI.x2.se, AVI.y2, code=3, angle=90, length=0.1, col="black")
points(AVI.x2, AVI.y2, pch=22, cex=1.4, col="black", bg="grey")


plot3<- plot(PAL.y3~PAL.x3, main="2009/10", cex.main=1.5, pch=21, cex=1.4, col="black", bg="black", xlim=c(-28.5,-22.5), ylim=c(7.0,9.2), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(AVI.y3~AVI.x3, pch=22, cex=1.4, col="grey", bg="grey")
points(CHA.y3~CHA.x3, pch=24, cex=1.4, col="slategrey", bg="slategrey")

arrows(PAL.x3, PAL.y3-PAL.y3.se, PAL.x3, PAL.y3+PAL.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(PAL.x3-PAL.x3.se, PAL.y3, PAL.x3+PAL.x3.se, PAL.y3, code=3, angle=90, length=0.1, col="black")
points(PAL.x3, PAL.y3, pch=21, cex=1.4, col="black", bg="black")

arrows(AVI.x3, AVI.y3-AVI.y3.se, AVI.x3, AVI.y3+AVI.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(AVI.x3-AVI.x3.se, AVI.y3, AVI.x3+AVI.x3.se, AVI.y3, code=3, 
       angle=90, length=0.1, col="black")
points(AVI.x3, AVI.y3, pch=22, cex=1.4, col="black", bg="grey")

arrows(CHA.x3, CHA.y3-CHA.y3.se, CHA.x3, CHA.y3+CHA.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(CHA.x3-CHA.x3.se, CHA.y3, CHA.x3+CHA.x3.se, CHA.y3, code=3, angle=90, length=0.1, col="black")
points(CHA.x3, CHA.y3, pch=24, cex=1.4, col="black", bg="slategrey")

dev.off()

#-----
# All stages graph with all years in one average

# Reset R, clear all objects
rm(list=ls()) 

#-----
# Load Data

#-----
setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

DataSet.ads.N<- read.csv("Ads.INS.DataSet.Sub3_FINAL.csv")

nrow(DataSet.ads.N)
DataSet.ads.N$Year.2<- as.factor(DataSet.ads.N$Year.2)
str(DataSet.ads.N)

range(DataSet.ads.N$delta.15.N)

#-----
list.files()

DataSet.ads.C<- read.csv("Ads.INS.DataSet.Sub4_FINAL.csv")

nrow(DataSet.ads.C)
DataSet.ads.C$Year.2<- as.factor(DataSet.ads.C$Year.2)
str(DataSet.ads.C)

range(DataSet.ads.C$delta.13.C)

#-----
list.files()

DataSet.cksD5.CN<- read.csv("Cks.D5.PAL.INS.DataSet.Sub2_FINAL.csv")

nrow(DataSet.cksD5.CN)
DataSet.cksD5.CN$Year.2<- as.factor(DataSet.cksD5.CN$Year.2)
str(DataSet.cksD5.CN)

range(DataSet.cksD5.CN$delta.15.N)
range(DataSet.cksD5.CN$delta.13.C)

#-----
list.files()

DataSet.cksD15.CN<- read.csv("Cks.D15.PAL.INS.DataSet.Sub2_FINAL.csv")

nrow(DataSet.cksD15.CN)
DataSet.cksD15.CN$Year.2<- as.factor(DataSet.cksD15.CN$Year.2)
str(DataSet.cksD15.CN)

range(DataSet.cksD15.CN$delta.15.N)
range(DataSet.cksD15.CN$delta.13.C)

#-----
list.files()

DataSet.cks5wk.PAL.CN<- read.csv("Cks.5wk.PAL.INS.DataSet.Sub7_FINAL.csv")

nrow(DataSet.cks5wk.PAL.CN)
DataSet.cks5wk.PAL.CN$Year.2<- as.factor(DataSet.cks5wk.PAL.CN$Year.2)
str(DataSet.cks5wk.PAL.CN)

range(DataSet.cks5wk.PAL.CN$delta.15.N)
range(DataSet.cks5wk.PAL.CN$delta.13.C)

#-----
list.files()

DataSet.cks5wk.WAP.N<- read.csv("Cks.5wk.WAP.INS.DataSet.Sub3_FINAL.csv")

nrow(DataSet.cks5wk.WAP.N)
DataSet.cks5wk.WAP.N$Year.2<- as.factor(DataSet.cks5wk.WAP.N$Year.2)
str(DataSet.cks5wk.WAP.N)

range(DataSet.cks5wk.WAP.N$delta.15.N)

#-----
list.files()

DataSet.cks5wk.WAP.C<- read.csv("Cks.5wk.WAP.INS.DataSet.Sub4_FINAL.csv")

nrow(DataSet.cks5wk.WAP.C)
DataSet.cks5wk.WAP.C$Year.2<- as.factor(DataSet.cks5wk.WAP.C$Year.2)
str(DataSet.cks5wk.WAP.C)

range(DataSet.cks5wk.WAP.C$delta.13.C)

#-----
# Subset DataSets into spp and stage, by isotope

#-----
# Adults

Ads.N.Spp1<- subset(DataSet.ads.N,Spp=="ADPE")
Ads.N.Spp2<- subset(DataSet.ads.N,Spp=="CHPE")
Ads.N.Spp3<- subset(DataSet.ads.N,Spp=="GEPE")

Ads.C.Spp1<- subset(DataSet.ads.C,Spp=="ADPE")
Ads.C.Spp2<- subset(DataSet.ads.C,Spp=="CHPE")
Ads.C.Spp3<- subset(DataSet.ads.C,Spp=="GEPE")

#-----
# Cks D5

CkD5.N.Spp1<- subset(DataSet.cksD5.CN,Spp=="ADPE")
CkD5.N.Spp2<- subset(DataSet.cksD5.CN,Spp=="CHPE")
CkD5.N.Spp3<- subset(DataSet.cksD5.CN,Spp=="GEPE")

CkD5.C.Spp1<- subset(DataSet.cksD5.CN,Spp=="ADPE")
CkD5.C.Spp2<- subset(DataSet.cksD5.CN,Spp=="CHPE")
CkD5.C.Spp3<- subset(DataSet.cksD5.CN,Spp=="GEPE")

#-----
# Cks D15

CkD15.N.Spp1<- subset(DataSet.cksD15.CN,Spp=="ADPE")
CkD15.N.Spp2<- subset(DataSet.cksD15.CN,Spp=="CHPE")
CkD15.N.Spp3<- subset(DataSet.cksD15.CN,Spp=="GEPE")

CkD15.C.Spp1<- subset(DataSet.cksD15.CN,Spp=="ADPE")
CkD15.C.Spp2<- subset(DataSet.cksD15.CN,Spp=="CHPE")
CkD15.C.Spp3<- subset(DataSet.cksD15.CN,Spp=="GEPE")

#-----
# Cks 5wk PAL

Ck5wkP.N.Spp1<- subset(DataSet.cks5wk.PAL.CN,Spp=="ADPE")
Ck5wkP.N.Spp2<- subset(DataSet.cks5wk.PAL.CN,Spp=="CHPE")
Ck5wkP.N.Spp3<- subset(DataSet.cks5wk.PAL.CN,Spp=="GEPE")

Ck5wkP.C.Spp1<- subset(DataSet.cks5wk.PAL.CN,Spp=="ADPE")
Ck5wkP.C.Spp2<- subset(DataSet.cks5wk.PAL.CN,Spp=="CHPE")
Ck5wkP.C.Spp3<- subset(DataSet.cks5wk.PAL.CN,Spp=="GEPE")

#-----
# Cks 5wk WAP

Ck5wkW.N.Loc1<- subset(DataSet.cks5wk.WAP.N,Location=="PAL")
Ck5wkW.N.Loc2<- subset(DataSet.cks5wk.WAP.N,Location=="AVI")
Ck5wkW.N.Loc3<- subset(DataSet.cks5wk.WAP.N,Location=="CHA")

Ck5wkW.C.Loc1<- subset(DataSet.cks5wk.WAP.C,Location=="PAL")
Ck5wkW.C.Loc2<- subset(DataSet.cks5wk.WAP.C,Location=="AVI")
Ck5wkW.C.Loc3<- subset(DataSet.cks5wk.WAP.C,Location=="CHA")

#-----
# Averages for each spp by stage for each isotope, and standard error. The variance estimate (var) is a built in function in R. The stdev is the sqrt of the variance. The se function was checked in excel and gave the same value.

variance<- function(x) sum((x-mean(x))^2)/(length(x)-1)
stdev<- function(x) sqrt(sum((x-mean(x))^2)/(length(x)-1))
se<- function(x) sqrt(var(x)/length(x))

# ADPE

AP.ads.mod1<- mean(Ads.C.Spp1$delta.13.C)
AP.ads.mod1.se<- se(Ads.C.Spp1$delta.13.C)

AP.ads.mod2<- mean(Ads.N.Spp1$delta.15.N)
AP.ads.mod2.se<- se(Ads.N.Spp1$delta.15.N)

AP.ckD5.mod3<- mean(CkD5.C.Spp1$delta.13.C)
AP.ckD5.mod3.se<- se(CkD5.C.Spp1$delta.13.C)

AP.ckD5.mod4<- mean(CkD5.N.Spp1$delta.15.N)
AP.ckD5.mod4.se<- se(CkD5.N.Spp1$delta.15.N)

AP.ckD15.mod5<- mean(CkD15.C.Spp1$delta.13.C)
AP.ckD15.mod5.se<- se(CkD15.C.Spp1$delta.13.C)

AP.ckD15.mod6<- mean(CkD15.N.Spp1$delta.15.N)
AP.ckD15.mod6.se<- se(CkD15.N.Spp1$delta.15.N)

AP.5wkP.mod7<- mean(Ck5wkP.C.Spp1$delta.13.C)
AP.5wkP.mod7.se<- se(Ck5wkP.C.Spp1$delta.13.C)

AP.5wkP.mod8<- mean(Ck5wkP.N.Spp1$delta.15.N)
AP.5wkP.mod8.se<- se(Ck5wkP.N.Spp1$delta.15.N)

AP.5wkW.mod9<- mean(Ck5wkW.C.Loc2$delta.13.C)
AP.5wkW.mod9.se<- se(Ck5wkW.C.Loc2$delta.13.C)

AP.5wkW.mod10<- mean(Ck5wkW.N.Loc2$delta.15.N)
AP.5wkW.mod10.se<- se(Ck5wkW.N.Loc2$delta.15.N)

AP.5wkW.mod11<- mean(Ck5wkW.C.Loc3$delta.13.C)
AP.5wkW.mod11.se<- se(Ck5wkW.C.Loc3$delta.13.C)

AP.5wkW.mod12<- mean(Ck5wkW.N.Loc3$delta.15.N)
AP.5wkW.mod12.se<- se(Ck5wkW.N.Loc3$delta.15.N)


#-----
# CHPE

CP.ads.mod1<- mean(Ads.C.Spp2$delta.13.C)
CP.ads.mod1.se<- se(Ads.C.Spp2$delta.13.C)

CP.ads.mod2<- mean(Ads.N.Spp2$delta.15.N)
CP.ads.mod2.se<- se(Ads.N.Spp2$delta.15.N)

CP.ckD5.mod3<- mean(CkD5.C.Spp2$delta.13.C)
CP.ckD5.mod3.se<- se(CkD5.C.Spp2$delta.13.C)

CP.ckD5.mod4<- mean(CkD5.N.Spp2$delta.15.N)
CP.ckD5.mod4.se<- se(CkD5.N.Spp2$delta.15.N)

CP.ckD15.mod5<- mean(CkD15.C.Spp2$delta.13.C)
CP.ckD15.mod5.se<- se(CkD15.C.Spp2$delta.13.C)

CP.ckD15.mod6<- mean(CkD15.N.Spp2$delta.15.N)
CP.ckD15.mod6.se<- se(CkD15.N.Spp2$delta.15.N)

CP.5wkP.mod7<- mean(Ck5wkP.C.Spp2$delta.13.C)
CP.5wkP.mod7.se<- se(Ck5wkP.C.Spp2$delta.13.C)

CP.5wkP.mod8<- mean(Ck5wkP.N.Spp2$delta.15.N)
CP.5wkP.mod8.se<- se(Ck5wkP.N.Spp2$delta.15.N)

#-----
# GEPE

GP.ads.mod1<- mean(Ads.C.Spp3$delta.13.C)
GP.ads.mod1.se<- se(Ads.C.Spp3$delta.13.C)

GP.ads.mod2<- mean(Ads.N.Spp3$delta.15.N)
GP.ads.mod2.se<- se(Ads.N.Spp3$delta.15.N)

GP.ckD5.mod3<- mean(CkD5.C.Spp3$delta.13.C)
GP.ckD5.mod3.se<- se(CkD5.C.Spp3$delta.13.C)

GP.ckD5.mod4<- mean(CkD5.N.Spp3$delta.15.N)
GP.ckD5.mod4.se<- se(CkD5.N.Spp3$delta.15.N)

GP.ckD15.mod5<- mean(CkD15.C.Spp3$delta.13.C)
GP.ckD15.mod5.se<- se(CkD15.C.Spp3$delta.13.C)

GP.ckD15.mod6<- mean(CkD15.N.Spp3$delta.15.N)
GP.ckD15.mod6.se<- se(CkD15.N.Spp3$delta.15.N)

GP.5wkP.mod7<- mean(Ck5wkP.C.Spp3$delta.13.C)
GP.5wkP.mod7.se<- se(Ck5wkP.C.Spp3$delta.13.C)

GP.5wkP.mod8<- mean(Ck5wkP.N.Spp3$delta.15.N)
GP.5wkP.mod8.se<- se(Ck5wkP.N.Spp3$delta.15.N)

#-----
# ADPE. Plot each stage by isotopes

AP.x1<- AP.ads.mod1
AP.y1<- AP.ads.mod2
AP.x1.se<- AP.ads.mod1.se
AP.y1.se<- AP.ads.mod2.se

AP.x2<- AP.ckD5.mod3
AP.y2<- AP.ckD5.mod4
AP.x2.se<- AP.ckD5.mod3.se
AP.y2.se<- AP.ckD5.mod4.se

AP.x3<- AP.ckD15.mod5
AP.y3<- AP.ckD15.mod6
AP.x3.se<- AP.ckD15.mod5.se
AP.y3.se<- AP.ckD15.mod5.se

AP.x4<- AP.5wkP.mod7
AP.y4<- AP.5wkP.mod8
AP.x4.se<- AP.5wkP.mod7.se
AP.y4.se<- AP.5wkP.mod8.se

AP.x5<- AP.5wkW.mod9
AP.y5<- AP.5wkW.mod10
AP.x5.se<- AP.5wkW.mod9.se
AP.y5.se<- AP.5wkW.mod10.se

AP.x6<- AP.5wkW.mod11
AP.y6<- AP.5wkW.mod12
AP.x6.se<- AP.5wkW.mod11.se
AP.y6.se<- AP.5wkW.mod12.se

#-----
# CHPE. Plot each year by isotopes

CP.x1<- CP.ads.mod1
CP.y1<- CP.ads.mod2
CP.x1.se<- CP.ads.mod1.se
CP.y1.se<- CP.ads.mod2.se

CP.x2<- CP.ckD5.mod3
CP.y2<- CP.ckD5.mod4
CP.x2.se<- CP.ckD5.mod3.se
CP.y2.se<- CP.ckD5.mod4.se

CP.x3<- CP.ckD15.mod5
CP.y3<- CP.ckD15.mod6
CP.x3.se<- CP.ckD15.mod5.se
CP.y3.se<- CP.ckD15.mod5.se

CP.x4<- CP.5wkP.mod7
CP.y4<- CP.5wkP.mod8
CP.x4.se<- CP.5wkP.mod7.se
CP.y4.se<- CP.5wkP.mod8.se

#-----
# GEPE. Plot each year by isotopes

GP.x1<- GP.ads.mod1
GP.y1<- GP.ads.mod2
GP.x1.se<- GP.ads.mod1.se
GP.y1.se<- GP.ads.mod2.se

GP.x2<- GP.ckD5.mod3
GP.y2<- GP.ckD5.mod4
GP.x2.se<- GP.ckD5.mod3.se
GP.y2.se<- GP.ckD5.mod4.se

GP.x3<- GP.ckD15.mod5
GP.y3<- GP.ckD15.mod6
GP.x3.se<- GP.ckD15.mod5.se
GP.y3.se<- GP.ckD15.mod5.se

GP.x4<- GP.5wkP.mod7
GP.y4<- GP.5wkP.mod8
GP.x4.se<- GP.5wkP.mod7.se
GP.y4.se<- GP.5wkP.mod8.se

#-----
# Plot

tiff(file="AllStage.CN.tiff",res=150,width=30,height=20,unit="cm")
par(mfrow=c(1,4))

plot1<- plot(AP.y1~AP.x1, main="Adults", cex.main=1.5, pch=16, col="black", xlim=c(-27.5,-22.5), ylim=c(7.25,9.5), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y1~CP.x1, pch=16, col="orangered")
points(GP.y1~GP.x1, pch=16, col="darkblue")

arrows(AP.x1, AP.y1-AP.y1.se, AP.x1, AP.y1+AP.y1.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x1-AP.x1.se, AP.y1, AP.x1+AP.x1.se, AP.y1, code=3, angle=90, length=0.1, col="black")
arrows(CP.x1, CP.y1-CP.y1.se, CP.x1, CP.y1+CP.y1.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x1-CP.x1.se, CP.y1, CP.x1+CP.x1.se, CP.y1, code=3, angle=90, length=0.1, col="orangered")
arrows(GP.x1, GP.y1-GP.y1.se, GP.x1, GP.y1+GP.y1.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x1-GP.x1.se, GP.y1, GP.x1+GP.x1.se, GP.y1, code=3, angle=90, length=0.1, col="darkblue")

plot2<- plot(AP.y2~AP.x2, main="Day 5 Chicks", cex.main=1.5, pch=16, col="black", xlim=c(-27.5,-22.5), ylim=c(7.25,9.5), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y2~CP.x2, pch=16, col="orangered")
points(GP.y2~GP.x2, pch=16, col="darkblue")

arrows(AP.x2, AP.y2-AP.y2.se, AP.x2, AP.y2+AP.y2.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x2-AP.x2.se, AP.y2, AP.x2+AP.x2.se, AP.y2, code=3, angle=90, length=0.1, col="black")
arrows(CP.x2, CP.y2-CP.y2.se, CP.x2, CP.y2+CP.y2.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x2-CP.x2.se, CP.y2, CP.x2+CP.x2.se, CP.y2, code=3, angle=90, length=0.1, col="orangered")
arrows(GP.x2, GP.y2-GP.y2.se, GP.x2, GP.y2+GP.y2.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x2-GP.x2.se, GP.y2, GP.x2+GP.x2.se, GP.y2, code=3, angle=90, length=0.1, col="darkblue")

plot3<- plot(AP.y3~AP.x3, main="Day 15 Chicks", cex.main=1.5, pch=16, col="black", xlim=c(-27.5,-22.5), ylim=c(7.25,9.5), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y3~CP.x3, pch=16, col="orangered")
points(GP.y3~GP.x3, pch=16, col="darkblue")

arrows(AP.x3, AP.y3-AP.y3.se, AP.x3, AP.y3+AP.y3.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x3-AP.x3.se, AP.y3, AP.x3+AP.x3.se, AP.y3, code=3, angle=90, length=0.1, col="black")
arrows(CP.x3, CP.y3-CP.y3.se, CP.x3, CP.y3+CP.y3.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x3-CP.x3.se, CP.y3, CP.x3+CP.x3.se, CP.y3, code=3, angle=90, length=0.1, col="orangered")
arrows(GP.x3, GP.y3-GP.y3.se, GP.x3, GP.y3+GP.y3.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x3-GP.x3.se, GP.y3, GP.x3+GP.x3.se, GP.y3, code=3, angle=90, length=0.1, col="darkblue")

plot4<- plot(AP.y4~AP.x4, main="Week 5 Chicks", cex.main=1.5, pch=16, col="black", xlim=c(-27.5,-22.5), ylim=c(7.25,9.5), xlab=expression(paste(delta^13,"C (‰)")), ylab=expression(paste(delta^15,"N (‰)")), cex.lab=1.75, cex.axis=1.2, par(mar=c(5.75,5.5,4.1,2.1)))
points(CP.y4~CP.x4, pch=16, col="orangered")
points(GP.y4~GP.x4, pch=16, col="darkblue")
points(AP.y5~AP.x5, pch=16, col="black")
points(AP.y6~AP.x6, pch=16, col="black")

arrows(AP.x4, AP.y4-AP.y4.se, AP.x4, AP.y4+AP.y4.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x4-AP.x4.se, AP.y4, AP.x4+AP.x4.se, AP.y4, code=3, angle=90, length=0.1, col="black")
arrows(CP.x4, CP.y4-CP.y4.se, CP.x4, CP.y4+CP.y4.se, code=3, angle=90, length=0.1, col="orangered")
arrows(CP.x4-CP.x4.se, CP.y4, CP.x4+CP.x4.se, CP.y4, code=3, angle=90, length=0.1, col="orangered")
arrows(GP.x4, GP.y4-GP.y4.se, GP.x4, GP.y4+GP.y4.se, code=3, angle=90, length=0.1, col="darkblue")
arrows(GP.x4-GP.x4.se, GP.y4, GP.x4+GP.x4.se, GP.y4, code=3, angle=90, length=0.1, col="darkblue")
arrows(AP.x5, AP.y5-AP.y5.se, AP.x5, AP.y5+AP.y5.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x5-AP.x5.se, AP.y5, AP.x5+AP.x5.se, AP.y5, code=3, angle=90, length=0.1, col="black")
arrows(AP.x6, AP.y6-AP.y6.se, AP.x6, AP.y6+AP.y6.se, code=3, angle=90, length=0.1, col="black")
arrows(AP.x6-AP.x6.se, AP.y6, AP.x6+AP.x6.se, AP.y6, code=3, angle=90, length=0.1, col="black")

dev.off()



