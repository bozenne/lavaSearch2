### sCorrect-getDesign.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:57) 
## Version: 
## Last-Updated: nov 19 2019 (10:52) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .getDesign
`.getDesign` <-
    function(object, data, ...) UseMethod(".getDesign")

## ** .getDesign.lm
.getDesign.lm <- function(object, data, add.Y, ...){
    X <- stats::model.matrix(stats::formula(object), droplevels(data)) ## drop unused factor (factor with 0 occurence)    
    attr(X,"assign") <- NULL
    attr(X,"contrasts") <- NULL
    if(add.Y){
        X <- cbind(X, data[endogenous(object, format = "long")])
    }
    return(X)
}

## ** .getDesign.gls
.getDesign.gls <- .getDesign.lm

## ** .getDesign.lme
.getDesign.lme <- .getDesign.lm

## ** .getDesign.lvmfit
.getDesign.lvmfit <- function(object, data, ...){
    return(as.matrix(data[,lava::manifest(object),drop=FALSE]))
}




######################################################################
### sCorrect-getDesign.R ends here
