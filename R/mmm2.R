### mmm2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:54) 
## Version: 
## Last-Updated: nov 29 2017 (17:31) 
##           By: Brice Ozenne
##     Update #: 37
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mmm2
### collect multiple marginal models
mmm2 <- function(...) {

     ret <- list(...)
     if (is.null(names(ret))) 
         names(ret) <- as.character(match.call(expand.dots = TRUE))[-1]
     class(ret) <- "mmm2"
     ret
}


### extract coefs
coef.mmm2 <- multcomp:::coef.mmm

## ### extract estimating functions
## estfun.mmm2 <- multcomp:::estfun.mmm

### set-up total covariance matrix
vcov.mmm2 <- function(object, return.null = TRUE,
                      adjust.residuals = TRUE, rescale = FALSE,
                      ...) {

    name.coef <- names(coef(object))
    n.coef <- length(name.coef)
    
    if(return.null){
        return(matrix(NA,n.coef,n.coef))
    }

    ## ** Extract influence functions from all models
    ls.iid <- lapply(object, iid2, adjust.residuals = adjust.residuals, power = 1/2)
    M.iid <- do.call(cbind, ls.iid)

    ## ** Rescale iid with sd/sd.robust 
    if(rescale){
        vec.sigma <- sqrt(unlist(lapply(object, function(x){diag(vcov(x))})))
        vec.sigma.robust <- sqrt(apply(M.iid^2,2,sum))
        M.iid <- sweep(M.iid, MARGIN = 2, FUN = "*", STATS = vec.sigma/vec.sigma.robust)
    }
    
    ## ** Compute covariance matrix
    robust.vcov <- crossprod(M.iid)

    ## ** export
    dimnames(robust.vcov) <- list(name.coef, name.coef)
    return(robust.vcov)
}

## * ls.lvmfit
vcov.ls.lvmfit <- vcov.mmm2

coef.ls.lvmfit <- coef.mmm2

##----------------------------------------------------------------------
### mmm2.R ends here
