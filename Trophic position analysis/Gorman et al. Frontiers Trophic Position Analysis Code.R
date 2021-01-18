#-----
# KB Gorman
# Pygoscelis Frontiers Manuscript
# Trophic Position Analysis (TPA)
# 17 Jan, 2020
#-----

# Clear workspace

rm(list=ls()) 

#-----
# Set Working Dir

setwd("/Users/kbgorman/Documents/KB Gorman Files/R Data Pengi")
getwd()
list.files()

#-----
# Load packages

install.packages("tRophicPosition")
library(tRophicPosition)
install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
install.packages("magrittr")
library(magrittr)

DataSet<- read.csv("PAL.5wk.ck.TL.Jan21-2_0708.csv")
DataSet<- read.csv("PAL.5wk.ck.TL.Jan21-2_0809.csv")
DataSet<- read.csv("PAL.5wk.ck.TL.Jan21-2_0910.csv")
DataSet<- read.csv("Cks.5wk.WAP.INS.DataSet.AVI_08.csv")
DataSet<- read.csv("Cks.5wk.WAP.INS.DataSet.AVI_09.csv")
DataSet<- read.csv("Cks.5wk.WAP.INS.DataSet.AVI_10.csv")
DataSet<- read.csv("Cks.5wk.WAP.INS.DataSet.CHA.csv")

DataSet
str(DataSet)

Penguin_plot<- screenFoodWeb(DataSet[which(DataSet$FG=="Penguin"),],
                             grouping =c("Spp", "FG"),
                             title = "Food web for Penguins (mean \u00B1 sd)",order = TRUE)

isotope_summary<- summariseIsotopeData(DataSet[which(DataSet$Location=="PAL"),],grouping =c("Spp", "FG"))
isotope_summary<- summariseIsotopeData(DataSet[which(DataSet$Location=="AVI"),],grouping =c("Spp", "FG"))
isotope_summary<- summariseIsotopeData(DataSet[which(DataSet$Location=="CHA"),],grouping =c("Spp", "FG"))
isotope_summary

write.table(isotope_summary, file="isotope_summary_PAL0708.csv", col.names=NA, sep=",")
write.table(isotope_summary, file="isotope_summary_PAL0809.csv", col.names=NA, sep=",")
write.table(isotope_summary, file="isotope_summary_PAL0910.csv", col.names=NA, sep=",")
write.table(isotope_summary, file="isotope_summary_AVI08.csv", col.names=NA, sep=",")
write.table(isotope_summary, file="isotope_summary_AVI09.csv", col.names=NA, sep=",")
write.table(isotope_summary, file="isotope_summary_AVI10.csv", col.names=NA, sep=",")
write.table(isotope_summary, file="isotope_summary_CHA10.csv", col.names=NA, sep=",")

DataSet_edited<- DataSet[DataSet$Location=="PAL",]
DataSet_edited<- DataSet[DataSet$Location=="AVI",]
DataSet_edited<- DataSet[DataSet$Location=="CHA",]

TDF_values <-TDF(author = "Post", element = "both")

PenguinList<- extractIsotopeData(DataSet_edited, 
                                 b1 = "Baseline",
                                 baselineColumn = "FG",
                                 consumersColumn = "Spp",
                                 groupsColumn = "Location",
                                 deltaC = TDF_values$deltaC,
                                 deltaN = TDF_values$deltaN)

summary(PenguinList)

cluster <- parallel::makePSOCKcluster(parallel::detectCores())

system.time(Penguin_TPmodels <- parallel::parLapply(cluster,
                                                        PenguinList,
                                                        lambda = 2.5,
                                                        multiModelTP,
                                                        adapt = 20000,
                                                        n.iter = 100000,
                                                        burnin = 20000,
                                                        n.chains = 5,
                                                        model = "oneBaseline"))

parallel::stopCluster(cluster)

model_summary<- fromParallelTP(Penguin_TPmodels, get = "summary")

write.table(model_summary, file="model_summary_PAL0708.csv", col.names=NA, sep=",")
write.table(model_summary, file="model_summary_PAL0809.csv", col.names=NA, sep=",")
write.table(model_summary, file="model_summary_PAL0910.csv", col.names=NA, sep=",")
write.table(model_summary, file="model_summary_AVI08.csv", col.names=NA, sep=",")
write.table(model_summary, file="model_summary_AVI09.csv", col.names=NA, sep=",")
write.table(model_summary, file="model_summary_AVI10.csv", col.names=NA, sep=",")
write.table(model_summary, file="model_summary_CHA10.csv", col.names=NA, sep=",")

ggplot_df <-fromParallelTP(Penguin_TPmodels, get = "summary")

ggplot_df2 <- dplyr::arrange(ggplot_df, group, mode)

values <-paste(ggplot_df2$group,"-",ggplot_df2$consumer)
ggplot_df2$species_ordered <-factor(values, levels = values, ordered = TRUE)

plot_penguins <-credibilityIntervals(ggplot_df2,
                                   x = "species_ordered",
                                   plotAlpha = FALSE,
                                   legend =c(0.92,0.15),
                                   group_by = "group",
                                   xlab = "Penguin species",
                                   scale_colour_manual =c("#2b83ba","#d7191c"), 
                                   print=TRUE)

plot_penguins <- plot_penguins+ggplot2::scale_x_discrete(labels =as.character(ggplot_df2$species))
print(plot_penguins)

TPs <-fromParallelTP(Penguin_TPmodels, get = "TP")

MTP <-list("PAL-ADPE" = TPs$`PAL-ADPE`,
           "PAL-CHPE" = TPs$`PAL-CHPE`,
           "PAL-GEPE" = TPs$`PAL-GEPE`)

MTP <-list("AVI-ADPE" = TPs$`AVI-ADPE`)

MTP <-list("CHA-ADPE" = TPs$`CHA-ADPE`)

(MTP_modes <-sapply(MTP, getPosteriorMode))

pairwiseComparisons(MTP, test = ">=")

pairwiseComparisons(MTP, test = "bhatt")


