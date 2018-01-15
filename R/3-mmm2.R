### mmm2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:54) 
## Version: 
## Last-Updated: jan 15 2018 (13:24) 
##           By: Brice Ozenne
##     Update #: 88
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


## * coef
#' @title Model Coefficients
#' @description Returns the model coefficients of a \code{mmm2} or \code{ls.lvmfit} object.
#' For internal use.
#' @name coef-multcomp
#'
#' @param object a \code{mmm2} or \code{ls.lvmfit} object.
#' @param ... Not used. Only for compatibility with the generic method.
#' 
#' @keywords internal

## * coef.mmm2
#' @rdname coef-multcomp
#' @method coef mmm2
#' @export
coef.mmm2 <- multcomp:::coef.mmm

## * coef.ls.lvmfit
#' @rdname coef-multcomp
#' @method coef ls.lvmfit
#' @export
coef.ls.lvmfit <- multcomp:::coef.mmm

## * vcov
#' @title Variance-Covariance Matrix for a Fitted Object
#' @description Returns the variance-covariance matrix of a fitted object.
#' For internal use.
#' @rdname vcov-multcomp
#' 
#' @param object a \code{mmm2} or \code{ls.lvmfit object}.
#' @param return.null if TRUE return a matrix filled with NA.
#' @param adjust.residuals small sample adjustement.
#' @param robust should robust standard error be used? Otherwise rescale the influence function with the standard error obtained from the information matrix.
#' @param ... Not used. Only for compatibility with the generic method.
#' @keywords internal


## * vcov.mm2
#' @rdname vcov-multcomp
#' @method vcov mmm2
#' 
#' @export
vcov.mmm2 <- function(object, return.null = TRUE,
                      adjust.residuals = TRUE, robust = TRUE,
                      ...) {

    name.coef <- names(stats::coef(object))
    n.coef <- length(name.coef)
    
    if(return.null){
        return(matrix(NA,n.coef,n.coef))
    }

    ## ** Extract influence functions from all models    
    ls.iid <- lapply(object, function(x){ ## x <- object[[1]]
        iIID <- iid2(x, adjust.residuals = adjust.residuals, power = 1/2)
        if(identical(class(x),"lm")){
            iIID <- iIID[,colnames(iIID)!="sigma2",drop=FALSE]                
        }
        return(iIID)
    })
    M.iid <- do.call(cbind, ls.iid)

    ## ** Rescale iid with sd/sd.robust 
    if(robust == FALSE){        
        ls.sd <- lapply(object, function(x){
            iVcov <- attr(residuals2(x, adjust.residuals = adjust.residuals,
                                     return.vcov.param = TRUE), "vcov.param")
            iSd <- sqrt(diag(iVcov))
            if(identical(class(x),"lm")){
                iSd <- iSd[names(iSd)!="sigma2"]                
            }
            return(iSd)
        })
        vec.sigma <- unlist(ls.sd)
        vec.sigma.robust <- sqrt(apply(M.iid^2,2,sum))
        M.iid <- sweep(M.iid, MARGIN = 2, FUN = "*", STATS = vec.sigma/vec.sigma.robust)
    }
    
    ## ** Compute covariance matrix
    robust.vcov <- crossprod(M.iid)

    ## ** export
    dimnames(robust.vcov) <- list(name.coef, name.coef)
    return(robust.vcov)
}

## * vcov.ls.lvmfit
#' @rdname vcov-multcomp
#' @method vcov ls.lvmfit
#' @export
vcov.ls.lvmfit <- vcov.mmm2


##----------------------------------------------------------------------
### mmm2.R ends here
