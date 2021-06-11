#----------------------------------------------------------------------------------------------------------------------------
#     R code 1) to perform random-effects Bayesian pairwise and network meta-analysis for aggregate continuous outcome data 
#                  <Normal likelihood, identity link, Random Effects> in Dias et al., 2013 in Appendix (PMID: 23104435)  
#                  Standardised Mean Difference after extending the aforementioned model
#                  One-stage pattern-mixture model with Informative Missingness Difference of Means (IMDoM) 
#                  (Mavridis et al., 2015 (PMID: 25393541)) under 25 different scenarios
#                  <Hierarchical, intervention-specific prior IMDoM> in Spineli et al., 2021 (PMID: 33406990)
#
#     R code 2) to perform random-effects Bayesian pairwise and network meta-analysis for aggregate binary outcome data 
#                  <Binomial likelihood, logit link, Random Effects> in Dias et al., 2013 in Appendix (PMID: 23104435)  
#                  One-stage pattern-mixture model with Informative Missingness Odds Ratio (IMOR) under 25 different scenarios 
#                  <Hierarchical, intervention-specific prior log IMOR> in Turner et al., 2015 (PMID: 25809313) & 
#                  Spineli et al., 2019 (PMID: 31018836)
#
#     Author: Loukia Spineli
#     Date: June 2021
#----------------------------------------------------------------------------------------------------------------------------



## Load libraries
list.of.packages <- c("dplyr", "R2jags")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("./R/data.preparation_function.R")
source("./R/sensitivity.analysis.mod_function.R")
source("./R/prepare.model_function.R")
source("./R/missingness.param.prior_function.R")
source("./R/heterogeneity.param.prior_function.R")



## Load datasets
load("./data/binary_PMA.RData") 
load("./data/binary_NMA.RData") 
load("./data/continuous_PMA.RData") 
load("./data/continuous_NMA.RData")



## Load empirical prior distributions for tau2
load("./data/tau2.priors.binary_NMA.RData") 
load("./data/tau2.priors.binary_PMA.RData")
load("./data/tau2.priors.continuous_PMA.RData") 
load("./data/tau2.priors.continuous_NMA.RData") 



## Indicate arguments for R2jags
D <- 0   # The direction of the outcome (i.e. harmful, or beneficial) does not matter 
n.chains <- 3
n.iter <- 100000
n.burnin <- 10000
n.thin <- 3




###########################################################################################################################
#
#                                  Run 108 PMAs on BINARY outcome data for each scenario                                                                                                                                                        
#
###########################################################################################################################


sens.binary.PMA <- list()
for(i in 1:length(unique(binary.PMA[, "ID"]))) {
  sens.binary.PMA[[i]] <- run.sensitivity(data = as.data.frame(binary.PMA[binary.PMA$ID == i, ]), 
                                          measure = "OR", 
                                          model = "RE", 
                                          assumption = "HIE-ARM", 
                                          heter.prior = list("lognormal", tau2.priors.binary.PMA[i, "mean"], 1/(tau2.priors.binary.PMA[i, "SD"])^2), 
                                          var.misspar = 1, 
                                          D = D, 
                                          n.chains = n.chains, 
                                          n.iter = n.iter, 
                                          n.burnin = n.burnin, 
                                          n.thin = n.thin)
}


## Save results on logOR in txt. format
LOR_PMA <- list()
for(i in 1:length(unique(binary.PMA[, "ID"]))){
  LOR_PMA[[i]] <- do.call(rbind, lapply(1:25, function(j) sens.binary.PMA[[i]]$EM))
}
(LOR.PMA <- do.call(rbind, lapply(1:length(unique(binary.PMA[, "ID"])), function(i) LOR_PMA[[i]])))
write.table(cbind(round(LOR.PMA, 4), rep(1:25, length(unique(binary.PMA[, "ID"]))), rep(1:length(unique(binary.PMA[, "ID"])), each = 25)), file = "./PMA_LOR.txt", sep = "\t", quote = F)




###########################################################################################################################
#
#                                   Run 29 NMAs on BINARY outcome data for each scenario                                                                                                                                                       
#
###########################################################################################################################


sens.binary.NMA <- item.bin.NMA <- list()
for(i in 1:length(binary.NMA)) {
  sens.binary.NMA[[i]] <- run.sensitivity(data = as.data.frame(binary.NMA[[i]]), 
                                          measure = "OR", 
                                          model = "RE", 
                                          assumption = "HIE-ARM", 
                                          heter.prior = list("lognormal", tau2.priors.binary.NMA[i, "Mean"], 1/(tau2.priors.binary.NMA[i, "SD"])^2), 
                                          var.misspar = 1, 
                                          D = D, 
                                          n.chains = n.chains, 
                                          n.iter = n.iter, 
                                          n.burnin = n.burnin, 
                                          n.thin = n.thin)
  item.bin.NMA[[i]] <- data.preparation(as.data.frame(binary.NMA[[i]]),
                                        measure = "OR")
}


## Save results on logOR in txt. format
LOR_NMA <- list()
for(i in 1:length(binary.NMA)){
  LOR_NMA[[i]] <- do.call(rbind, lapply(1:25, function(j) sens.binary.NMA[[i]]$EM))
}
(LOR.NMA <- do.call(rbind, lapply(1:length(binary.NMA), function(i) LOR_NMA[[i]])))
write.table(cbind(round(LOR.NMA, 4), unlist(lapply(1:length(binary.NMA), function(i) rep(1:25, each = item.bin.NMA[[i]]$nt*(item.bin.NMA[[i]]$nt - 1)*0.5))), unlist(lapply(1:length(binary.NMA), function(i) rep(i, each = item.bin.NMA[[i]]$nt*(item.bin.NMA[[i]]$nt - 1)*0.5*25)))), file = "./NMA_LOR.txt", sep = "\t", quote = F)




###########################################################################################################################
#
#                                 Run 13 PMAs on CONTINUOUS outcome data for each scenario                                                                                                                                                                                                                         
#
###########################################################################################################################

 
sens.continuous.PMA <- list()
for(i in 1:length(unique(continuous.PMA[, "ID"]))) {
  sens.continuous.PMA[[i]] <- run.sensitivity(data = as.data.frame(continuous.PMA[continuous.PMA$ID == i, ]), 
                                              measure = "SMD", 
                                              model = "RE", 
                                              assumption = "HIE-ARM", 
                                              heter.prior = list("logt", tau2.priors.continuous.PMA[i, "mean"], 1/(tau2.priors.continuous.PMA[i, "SD"])^2), 
                                              var.misspar = 1, 
                                              D = D, 
                                              n.chains = n.chains, 
                                              n.iter = n.iter, 
                                              n.burnin = n.burnin, 
                                              n.thin = n.thin)
}


## Save results on SMD in txt. format
SMD_PMA <- list()
for(i in 1:length(unique(continuous.PMA[, "ID"]))){
  SMD_PMA[[i]] <- do.call(rbind, lapply(1:25, function(j) sens.continuous.PMA[[i]]$EM))
}
(SMD.PMA <- do.call(rbind, lapply(1:length(unique(continuous.PMA[, "ID"])), function(i) SMD_PMA[[i]])))
write.table(cbind(round(SMD.PMA, 4), rep(1:25, length(unique(continuous.PMA[, "ID"]))), rep(1:length(unique(continuous.PMA[, "ID"])), each = 25)), file = "./PMA_SMD.txt", sep = "\t", quote = F)




###########################################################################################################################
#
#                             Run 5 NMAs on CONTINUOUS for each scenario (25 scenarios in total)                                                                                       
#
###########################################################################################################################


sens.continuous.NMA <- item.con.NMA <- list()
for(i in 1:length(continuous.NMA)) {
  sens.continuous.NMA[[i]] <- run.sensitivity(data = as.data.frame(continuous.NMA[[i]]), 
                                              measure = "SMD", 
                                              model = "RE", 
                                              assumption = "HIE-ARM", 
                                              heter.prior = list("logt", tau2.priors.continuous.NMA[i, "mean"], 1/(tau2.priors.continuous.NMA[i, "SD"])^2), 
                                              var.misspar = 1, 
                                              D = D, 
                                              n.chains = n.chains, 
                                              n.iter = n.iter, 
                                              n.burnin = n.burnin, 
                                              n.thin = n.thin)
  item.bin.NMA[[i]] <- data.preparation(as.data.frame(continuous.NMA[[i]]),
                                        measure = "SMD")
}


## Save results on SMD in txt. format
SMD_NMA <- list()
for(i in 1:length(continuous.NMA)){
  SMD_NMA[[i]] <- do.call(rbind, lapply(1:25, function(j) sens.continuous.NMA[[i]]$EM))
}
(SMD.NMA <- do.call(rbind, lapply(1:length(continuous.NMA), function(i) SMD_NMA[[i]])))
write.table(cbind(round(SMD.NMA, 4), unlist(lapply(1:length(continuous.NMA), function(i) rep(1:25, each = item.con.NMA[[i]]$nt*(item.con.NMA[[i]]$nt - 1)*0.5))), unlist(lapply(1:length(continuous.NMA), function(i) rep(i, each = item.con.NMA[[i]]$nt*(item.con.NMA[[i]]$nt - 1)*0.5*25)))), file = "./NMA_SMD.txt", sep = "\t", quote = F)



