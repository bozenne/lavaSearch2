### sCorrect-getIndexOmega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 25 2019 (10:52) 
## Version: 
## Last-Updated: feb 20 2020 (10:46) 
##           By: Brice Ozenne
##     Update #: 91
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Identify the Endogenous Variables
#' @description Identify the endogenous variables, i.e., returns a vector with length the number of observations,
#' whose values are the index of the repetitions.
#' @name getIndexOmega-internal
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param data dataset.
#' @param ... [internal] Only used by the generic method.
#'  
#' @concept extractor
#' @keywords internal
`.getIndexOmega` <-
    function(object, data, ...) UseMethod(".getIndexOmega")

## * Examples
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### linear model ####
#' e.lm <- lm(Y1~X1, data = dW)
#' .getIndexOmega(e.lm, data = dW)
#'
#' #### gls model ####
#' e.gls1 <- gls(Y1~X1, data = dW)
#' .getIndexOmega(e.gls1, data = dW)
#' 
#' e.gls2 <- gls(Y~X1, correlation = corCompSymm(form=~1|id), data = dL)
#' .getIndexOmega(e.gls2, data = dL)
#'
#' e.gls3 <- gls(Y~X1, correlation = corCompSymm(form=~time|id), data = dL)
#' .getIndexOmega(e.gls3, data = dL)
#'
#' e.gls4 <- gls(Y~X1, weight = varIdent(form=~1|time2), data = dL)
#' .getIndexOmega(e.gls4, data = dL)
#'
#' e.gls5 <- gls(Y~X1, weight = varIdent(form=~1|time2),
#'               correlation = corSymm(form=~time|id), data = dL)
#' .getIndexOmega(e.gls5, data = dL)
#'
#' #### lme model ####
#' e.lme <- lme(Y~X1, random=~1|id, data = dL)
#' .getIndexOmega(e.lme, data = dL)
#' 
#' #### lvm model ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' .getIndexOmega(e.lvm, data = dW)
#' 

## * .getIndexOmega.lm
#' @rdname getIndexOmega-internal
.getIndexOmega.lm <- function(object, data, ...){

    ## ** check missing value in X
    name.manifest <- lava::manifest(object)
    test.na <- Reduce("+",lapply(name.manifest, function(iCol){is.na(data[[iCol]])}))
    if(any(test.na>0)){
        stop("Does not support missing values in \"lm\" objects. \n",
             "Remove rows in the dataset with missing values. \n")
    }

    ## ** export
    n.obs <- NROW(data)
    return(rep(1,n.obs))
}

## * .getIndexOmega.gls
#' @rdname getIndexOmega-internal
.getIndexOmega.gls <- function(object, data, ...){

    ## ** check missing value in X or grouping variables
    name.manifest <- manifest(object, format = "long")
    test.na <- Reduce("+",lapply(name.manifest, function(iCol){is.na(data[[iCol]])}))
    if(any(test.na>0)){
        stop("Does not support missing values in \"lm\" objects. \n",
             "Remove rows in the dataset with missing values. \n")
    }

    ## ** ref.group
    if(!is.null(object$modelStruct$varStruct)){
        var.var <- all.vars(formula(object$modelStruct$varStruct))
        if(length(var.var)>1){
            stop("Can only handle one covariate in the formula for the variance (argument \'weight\') \n")
        }
        level.var.var <- attr(unclass(object$modelStruct$varStruct),"groupNames")
    }
    ## ** index.Omega
    n.obs <- NROW(data)
    if(!is.null(object$modelStruct$corStruct)){
        iFormula <- formula(object$modelStruct$corStruct)
        varIndex.cor <- all.vars(nlme::getCovariateFormula(iFormula))
        var.cor <- setdiff(all.vars(iFormula),varIndex.cor)
    }else{
        varIndex.cor <- NULL
        var.cor <- NULL
    }
    if(!is.null(object$modelStruct$reStruct)){
        iFormula <- formula(object$modelStruct$reStruct)
        varIndex.re <- sapply(iFormula,function(iF){all.vars(nlme::getCovariateFormula(iF))})
        var.re <- sapply(iFormula,function(iF){setdiff(all.vars(iF),all.vars(nlme::getCovariateFormula(iF)))})
    }else{
        varIndex.re <- NULL
        var.re <- NULL
    }
    
    if(!is.null(object$modelStruct$corStruct) && length(varIndex.cor)>0){
        index.Omega <- as.numeric(as.factor(data[[varIndex.cor]]))
        if(!is.null(object$modelStruct$varStruct)){
            index.Omega2 <- as.numeric(factor(data[[var.var]], levels = level.var.var))
            attr(index.Omega,"index2endogenous") <- tapply(index.Omega, index.Omega2, unique)
        }
    }else if(!is.null(object$modelStruct$varStruct)){
        index.Omega <- as.numeric(factor(data[[var.var]], levels = level.var.var))
    }else if(!is.null(object$modelStruct$reStruct)){
        id <- factor(data[[var.re]])
        index.Omega <- Reduce("+",lapply(levels(id), function(iL){
            (iL==id)*cumsum(iL==id)
        }))
    }else if(!is.null(object$modelStruct$corStruct)){
        id <- factor(data[[var.cor]])
        index.Omega <- Reduce("+",lapply(levels(id), function(iL){
            (iL==id)*cumsum(iL==id)
        }))
    }else{
        index.Omega <- rep(1,n.obs)
    }
    
    ## ** ref.group
    if(!is.null(object$modelStruct$varStruct)){
        attr(index.Omega,"ref") <- level.var.var[1]
    }
    if(is.null(attr(index.Omega,"index2endogenous"))){
        Uindex.Omega <- unique(index.Omega)
        attr(index.Omega,"index2endogenous") <- setNames(as.list(Uindex.Omega),Uindex.Omega)
    }
    return(index.Omega)
}


## * .getIndexOmega.lme
#' @rdname getIndexOmega-internal
.getIndexOmega.lme <- .getIndexOmega.gls

## * .getIndexOmega.lvm
#' @rdname getIndexOmega-internal
.getIndexOmega.lvm <- function(object, data, ...){

    ## ** check missing value in exogenous variables
    name.exogenous <- exogenous(object)
    missing.var <- name.exogenous[name.exogenous %in% names(data) == FALSE]

    if(length(missing.var)>0){
        cat2bin <- var2dummy(object$model, var = names(dW), data = dW)
        name.exogenous[name.exogenous %in% missing.var] <- names(cat2bin)[cat2bin %in% missing.var]
        name.exogenous <- unique(name.exogenous)
    }
    test.na <- Reduce("+",lapply(name.exogenous, function(iCol){is.na(data[[iCol]])}))
    if(any(test.na>0)){
        stop("Does not support missing values in \"lm\" objects. \n",
             "Remove rows in the dataset with missing values. \n")
    }

    ## ** index.Omega
    n.obs <- NROW(data)
    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)

    M.index <- matrix(1:n.endogenous, nrow = n.obs, ncol = n.endogenous, byrow = TRUE)
    index.na <- which(is.na(as.data.frame(data)[,name.endogenous]))
    if(length(index.na)>0){
        M.index[index.na] <- NA
    }
    return(as.double(t(M.index)))
}

## * .getIndexOmega.lvmfit
#' @rdname getIndexOmega-internal
.getIndexOmega.lvmfit <- .getIndexOmega.lvm

######################################################################
### sCorrect-getIndexOmega.R ends here
