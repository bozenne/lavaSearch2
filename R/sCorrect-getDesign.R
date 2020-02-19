### sCorrect-getDesign.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:57) 
## Version: 
## Last-Updated: feb 19 2020 (12:10) 
##           By: Brice Ozenne
##     Update #: 18
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

## ** .getDesign.lvm
.getDesign.lvm <- function(object, data, ...){
    if(all(apply(data[,lava::manifest(object),drop=FALSE],2,is.numeric))){
        return(as.matrix(data[,lava::manifest(object),drop=FALSE]))
    }else{
        stop("Cannot export design matrix because some of the variables are categorical. \n")
    }
}

## ** .getDesign.lvmfit
.getDesign.lvmfit <- function(object, data, ...){
    return(as.matrix(data[,lava::manifest(object),drop=FALSE]))
}



######################################################################
### sCorrect-getDesign.R ends here
