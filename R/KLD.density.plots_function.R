#----------------------------------------------------------------------------------------------------------------------------
#     A function to obtain the Kullback-Leibler divergence (KLD, comparing two normal distributions)
#     A function to create a panel of density plots on the effect measure and KLD under all missingness scenarios
#     Author: Loukia Spineli
#     Date: June 2021
#----------------------------------------------------------------------------------------------------------------------------



## Function for Kullback-Leibler Divergence (comparing two univariate normal distributions)
KLD.measure.univ <- function(mean.y, sd.y, mean.x, sd.x){
  
  # x is the 'truth' (e.g. the MAR assumption)
  KLD.xy <- 0.5*(((sd.x/sd.y)^2) + ((mean.y - mean.x)^2)/(sd.y^2) - 1 + 2*log(sd.y/sd.x))
  
  return(list(KLD.xy))
}



KLD.plots <- function(ES.mat, primary.scenar, compar, outcome, drug.names){

  
  if(outcome == "binary"){
    
    scenarios <- c("1/3", "1/2", "1", "2", "3")
    
  } else {
    
    scenarios <- c("-2", "-1", "0", "1", "2")
    
  }
  n.scenar <- length(scenarios)
  names(scenarios) <- as.character(1:n.scenar)
  

  (nt <- (1 + sqrt(1 + 8*(length(ES.mat[, 1])/n.scenar^2)))/2)  # The quadratic formula for the roots of the general quadratic equation
  ## Group the results per scenario
  upper.ES <- lower.ES <- sd.ES <- ES <- matrix(NA, nrow = n.scenar^2, ncol = (nt*(nt - 1))/2)  
  for(i in 1:(n.scenar^2)){
    for(j in 1:(nt*(nt - 1))/2){
      ES[i, j] <- ES.mat[j + (nt*(nt - 1)/2)*(i - 1), 1]
      sd.ES[i, j] <- ES.mat[j + (nt*(nt - 1)/2)*(i - 1), 2]
      lower.ES[i, j] <- ES.mat[j + (nt*(nt - 1)/2)*(i - 1), 3]
      upper.ES[i, j] <- ES.mat[j + (nt*(nt - 1)/2)*(i - 1), 4]
    }
  }
  
  comparison <- matrix(combn(drug.names, 2), nrow = length(combn(drug.names, 2))/2, ncol = 2, byrow = T)
  
  KLD.xy <- list()
  
  for(i in 1:(nt*(nt - 1)/2)){
    
    KLD.xy[[i]] <- rep(NA, n.scenar*n.scenar)
    
    for(j in 1:(n.scenar*n.scenar)[-primary.scenar]){
      
      ## Returns the KLD of informative scenario j when compared with MAR (MAR as 'true') for reference-comparison i  
      KLD.xy[[i]][j] <- KLD.measure.univ(ES[j, i], sd.ES[j, i], ES[primary.scenar, i], sd.ES[primary.scenar, i])[[1]]
      
    }
    
    KLD.xy[[i]][primary.scenar] <- 0
  }
  
  mean.x <- rep(ES[primary.scenar, compar], n.scenar*n.scenar); sd.x <- rep(sd.ES[primary.scenar, compar], n.scenar*n.scenar)
  
  time0 <- prob <- prob.y0 <- scen.all <- IL0 <- KLD0 <- lower0 <- upper0 <- list()  
  for(i in 1:(n.scenar*n.scenar)){
    time0[[i]] <- seq(ES[i, compar] - 3.3*sd.ES[i, compar], ES[i, compar] + 3.3*sd.ES[i, compar], 0.1)
    prob[[i]] <- dnorm(time0[[i]], ES[i, compar], sd.ES[i, compar]) - dnorm(time0[[i]], mean.x[i], sd.x[i])
    prob.y0[[i]] <- dnorm(time0[[i]], ES[i, compar], sd.ES[i, compar])
    scen.all[[i]] <- rep(cbind(rep(1:5, each = 5), rep(1:5, 5))[i, ], length(time0[[i]])) # first column 'active', second column 'control'
    KLD0[[i]] <- rep(round(KLD.xy[[compar]][i], 3), length(time0[[i]]))
    lower0[[i]] <- rep(lower.ES[i, compar], length(time0[[i]]))
    upper0[[i]] <- rep(upper.ES[i, compar], length(time0[[i]]))
  }
  out <- as.data.frame(cbind(unlist(time0), unlist(prob), unlist(prob.y0), matrix(unlist(scen.all), nrow = length(unlist(scen.all))/2, ncol = 2, byrow = T), 
                             unlist(KLD0), unlist(lower0), unlist(upper0)))
  colnames(out) <- c("time", "x", "prob.y", "exper", "ctrl", "KLD", "lower", "upper")
  
  p <- ggplot(data = out, aes(time, x)) + 
         geom_rect(aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf), fill = "grey85") +
         stat_function(fun = function(z){dnorm(z, mean.x, sd.x)}, col = "#D55E00", size = 1.3) +
         geom_line(aes(time, x), size = 1.5, col = "blue") +
         geom_line(aes(time, prob.y), size = 1.5, col = "black") +
         geom_area(fill = "lightblue") +
         geom_vline(xintercept = 0, linetype = 2) + 
         geom_hline(yintercept = 0) +
         geom_text(x = -Inf, y = Inf, aes(label = paste0("KLD=", KLD)), size = 3.5, hjust = -0.2, vjust = 1.4, check_overlap = T) +
         facet_grid(ctrl ~ exper, switch = "y", labeller = labeller(ctrl = scenarios, exper = scenarios)) +
         theme_test() +
         labs(x = ifelse(outcome == "binary", paste("log odds ratio of", comparison[compar, 2], "versus", comparison[compar, 1]),
                         paste("standardised mean difference of", comparison[compar, 2], "versus", comparison[compar, 1])), y = "") +
         theme(axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 13), axis.text.y = element_text(size = 13),
               strip.text = element_text(size = 13), strip.background = element_rect(color = "grey30", fill = "grey90"),
               strip.placement = "outside")
  return(p)
}
