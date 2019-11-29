### manifest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 19 2019 (11:06) 
## Version: 
## Last-Updated: nov 25 2019 (11:54) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

manifest.lm <- function(x,...){
    return(c(endogenous(x),exogenous(x)))
}

manifest.gls <- function(x,...){
    out <- c(endogenous(x, ...),exogenous(x))
    if(!is.null(x$modelStruct$reStruct)){
        out <- union(out,sapply(formula(x$modelStruct$reStruct), all.vars))
    }
    if(!is.null(x$modelStruct$corStruct)){
        out <- union(out,all.vars(formula(x$modelStruct$corStruct)))
    }
    if(!is.null(x$modelStruct$varStruct)){
        out <- union(out,all.vars(formula(x$modelStruct$varStruct)))
    }

    return(out)
}

manifest.lme <- manifest.gls

######################################################################
### manifest.R ends here
