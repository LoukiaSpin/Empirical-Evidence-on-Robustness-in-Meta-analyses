#----------------------------------------------------------------------------------------------------------------------------
#     A function to determine the 'mode' of the mean value the informative missingness parameter 
#     Author: Loukia Spineli
#     Date: June 2021
#----------------------------------------------------------------------------------------------------------------------------


missingness.param.prior <- function(assumption, mean.misspar) {


  ## Condition regarding the specification of the prior mean ('mean.misspar') for the missingness parameter
  if (missing(mean.misspar) & (is.element(assumption, c("HIE-ARM", "IDE-ARM" )))) {

    mean.misspar <- rep(0.0001, 2)

  }  else if (missing(mean.misspar) & (!is.element(assumption, c("HIE-ARM", "IDE-ARM")))) {

    mean.misspar <- 0.0001

  }  else if (!missing(mean.misspar) & (is.element(assumption, c("HIE-ARM", "IDE-ARM"))) & is.null(dim(mean.misspar))) {

    mean.misspar <- rep(ifelse(mean.misspar == 0, 0.0001, mean.misspar), 2)

  } else if (!missing(mean.misspar) & (is.element(assumption, c("HIE-ARM", "IDE-ARM"))) & !is.null(dim(mean.misspar))) {

    mean.misspar <- as.vector(mean.misspar)
    mean.misspar[1] <- ifelse(mean.misspar[1] == 0, 0.0001, mean.misspar[1])
    mean.misspar[2] <- ifelse(mean.misspar[2] == 0, 0.0001, mean.misspar[2])

  } else {

    mean.misspar <- mean.misspar

  }

}
