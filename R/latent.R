### latent.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (13:48) 
## Version: 
## Last-Updated: dec  9 2019 (10:44) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

latent.lm <- function (x, ...){
    return(NULL)
}
latent.gls <- latent.lm
latent.lme <- function(x, ...){
    return(paste0("eta.",names(ranef(x))))
}

######################################################################
### latent.R ends here
