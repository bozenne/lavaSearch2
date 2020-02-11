### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  3 2020 (18:29) 
## Version: 
## Last-Updated: feb 11 2020 (17:35) 
##           By: Brice Ozenne
##     Update #: 15
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
        dtransform <- function(x){numDeriv::jacobian(transform, x)[1,1]}
    }
    object[,"std.error"] <- object[,"std.error"]*dtransform(object[,"estimate"])
    object[,"estimate"] <- transform(object[,"estimate"])
    if("ci.lower" %in% names(object)){
        object[,"ci.lower"] <- transform(object[,"ci.lower"])
    }
    if("ci.upper" %in% names(object)){
        object[,"ci.upper"] <- transform(object[,"ci.upper"])
    }
    return(object)
}

######################################################################
### transform.R ends here
