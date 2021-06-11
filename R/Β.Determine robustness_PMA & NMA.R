#----------------------------------------------------------------------------------------------------------------------------
#     R code 1) to calculate the robustness index for all analyses from 'Part A'
#            2) create the heatmap of robustness index for each possible comparison of the network, and 
#            3) create the density plots for the primary and alternative re-analysis, and the Kullback-Leibler Divergence
#
#     Author: Loukia Spineli
#     Date: June 2021
#----------------------------------------------------------------------------------------------------------------------------



## Load libraries
list.of.packages <- c("ggplot2", "reshape2")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("./R/robustness.index_function.R")
source("./R/heatmap.robustness_function.R")
source("./R/KLD.density.plots_function.R")
source("./R/possible.comparisons.id_function.R")



## Load JAGS results on logOR and SMDs (PMA & NMA)
bin.pma <- read.table("./data/PMA_LOR.txt", header = T)
con.pma <- read.table("./data/PMA_SMD.txt", header = T)
bin.nma <- read.table("./data/NMA_LOR.txt", header = T) 
con.nma <- read.table("./data/NMA_SMD.txt", header = T)




#########################################################################################################################
# 
#                                     Calculate the Robustness Index for each PMA                                       
#
#########################################################################################################################


## Binary outcomes
RI.bin.pma <- rep(NA, length(unique(bin.pma$MA)))
for (i in 1:length(unique(bin.pma$MA))) {
  RI.bin.pma[i] <- RobustnessIndex(bin.pma[bin.pma$MA == i, 2:3], 
                                   threshold = 0.28, 
                                   primary.scenar = 13,
                                   nt = 2)$robust
}
table(RI.bin.pma[-81])    # Exclude PMA (on binary outcome) no 81 due to non-convergence!


## Continuous outcomes
RI.con.pma <- rep(NA, length(unique(con.pma$MA)))
for (i in 1:length(unique(con.pma$MA))) {
  RI.con.pma[i] <- RobustnessIndex(con.pma[con.pma$MA == i, 2:3], 
                                   threshold = 0.17, 
                                   primary.scenar = 13,
                                   nt = 2)$robust
}
table(RI.con.pma)




#########################################################################################################################
# 
#                                     Calculate the Robustness Index for each NMA                                       
#
#########################################################################################################################


## Binary outcomes
RI.bin.nma0 <- list()
RI.bin.nma <- rep(NA, length(unique(bin.nma$NMA)))
for (i in 1:length(unique(bin.nma$NMA))) {
  RI.bin.nma0[[i]] <- RobustnessIndex(bin.nma[bin.nma$NMA == i, 2:3], 
                                      threshold = 0.28, 
                                      primary.scenar = 13, 
                                      nt = (1 + sqrt(1 + 8*(length(bin.nma[bin.nma$NMA == i, 1])/25)))/2)$robust
  RI.bin.nma[i] <- ifelse(length(unique(RI.bin.nma0[[i]])) > 1, "frail", unique(RI.bin.nma0[[i]]))
}
table(RI.bin.nma[-22])    # Exclude NMA (on binary outcome) no 22 due to non-convergence!


## Continuous outcomes
RI.con.nma0 <- list()
RI.con.nma <- rep(NA, length(unique(con.nma$NMA)))
for (i in 1:length(unique(con.nma$NMA))) {
  RI.con.nma0[[i]] <- RobustnessIndex(con.nma[con.nma$NMA == i, 2:3], 
                                      threshold = 0.17, 
                                      primary.scenar = 13, 
                                      nt = (1 + sqrt(1 + 8*(length(con.nma[con.nma$NMA == i, 1])/25)))/2)$robust
  RI.con.nma[i] <- ifelse(length(unique(RI.con.nma0[[i]])) > 1, "frail", unique(RI.con.nma0[[i]]))
}
table(RI.con.nma)




#######################################################################################################################
#
#                        Create the heatmap of robustness index and the panel of density plots                                              
#                                using the network of Liu et al., 2013 (PMID: 24098546)                                 
#
#######################################################################################################################


## Obtain the robustness index for each possible comparison of the network
RI <- RobustnessIndex(bin.nma[bin.nma$NMA == 14, 2:3], 
                threshold = 0.28, 
                primary.scenar = 13, 
                nt = (1 + sqrt(1 + 8*(length(bin.nma[bin.nma$NMA == 14, 1])/25)))/2)$RI


## Heatmap of robustness index for the network
HeatMap.AllComparisons.RI(RI, 
                          drug.names = c("placebo", "pramipexole", "SNRI", "SSRI", "TCA", "pergolide"), 
                          threshold = 0.28)


## Find the 'id' of the comparison of interest to be used in 'compar' of the 'KLD.plots' function
possible.comparisons.id(drug.names = c("placebo", "pramipexole", "SNRI", "SSRI", "TCA", "pergolide"))


## Panel of density plots on the effect measure under all scenarios and illustration of the Kullback-Leibler Divergence
KLD.plots(bin.nma[bin.nma$NMA == 14, 2:5], 
          primary.scenar = 13, 
          compar = 3,           # This refers to the comparison of SNRI with placebo
          outcome = "binary", 
          drug.names = c("placebo", "pramipexole", "SNRI", "SSRI", "TCA", "pergolide"))



