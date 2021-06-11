#----------------------------------------------------------------------------------------------------------------------------
#     A function to identify the id-number of the pairwise comparison of interest 
#     Author: Loukia Spineli
#     Date: June 2021
#----------------------------------------------------------------------------------------------------------------------------


possible.comparisons.id <- function(drug.names) {


  ## Obtain all unique pairwise comparisons using the 'combn' functions
  nt <- length(drug.names)
  poss.pair.comp <- data.frame(t(combn(1:nt, 2))[, 2], t(combn(1:nt, 2))[, 1])
  poss.pair.comp$comp <- paste0(poss.pair.comp[, 1], "vs", poss.pair.comp[, 2])
  colnames(poss.pair.comp) <- c("treat1", "treat2", "comp")


  poss.pair.comp.name <- poss.pair.comp
  ## Replace intervention id with their original name
  # All possible comparisons - Treat1 (non-baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp[, 1])))) {
    poss.pair.comp.name[poss.pair.comp.name$treat1 == i, 1] <- drug.names[i]
  }
  # Observed comparisons - Treat2 (baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp[, 2])))) {
    poss.pair.comp.name[poss.pair.comp.name$treat2 == i, 2] <- drug.names[i]
  }


  return(list(poss.comp = data.frame(ID = 1:length(poss.pair.comp$comp), poss.pair.comp[, 1:2], poss.pair.comp.name[, 1:2])))

}
