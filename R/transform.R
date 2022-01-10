### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  3 2020 (18:29) 
## Version: 
## Last-Updated: Jan  5 2022 (12:47) 
##           By: Brice Ozenne
##     Update #: 20
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @export
transformSummaryTable <- function(object, transform = NULL, conf.level = 0.95){
    if(is.null(transform)){
        return(object)
    }else if(identical(transform,"atanh")){
        transform <- atanh
        dtransform <- function(x){1/(1-x^2)}
    }else if(identical(transform,"exp")){
        transform <- exp
        dtransform <- function(x){exp(x)}
    }else if(identical(transform,"log")){
        transform <- log
        dtransform <- function(x){1/x}
    }else if(identical(transform,"loglog")){
        transform <- function(x){log(log(x))}
        dtransform <- function(x){1/(-x*log(x))}
    }else if(identical(transform,"cloglog")){
        transform <- function(x){log(log(1-x))}
        dtransform <- function(x){1/(-(1-x)*log(1-x))}
    }else if(!is.null(attr(transform,"derivative"))){
        dtransform <- attr(transform,"derivative")
    }else{
        dtransform <- function(x){diag(numDeriv::jacobian(transform, x))}
    }
    object[,"se"] <- object[,"se"]*dtransform(object[,"estimate"])
    object[,"estimate"] <- transform(object[,"estimate"])
    if("lower" %in% names(object)){
        object[,"lower"] <- transform(object[,"lower"])
    }
    if("upper" %in% names(object)){
        object[,"upper"] <- transform(object[,"upper"])
    }
    return(object)
}

######################################################################
### transform.R ends here
