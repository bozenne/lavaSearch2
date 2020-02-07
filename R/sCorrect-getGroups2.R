### sCorrect-getGroups2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:58) 
## Version: 
## Last-Updated: feb  6 2020 (17:12) 
##           By: Brice Ozenne
##     Update #: 141
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Reconstruct the Cluster variable
#' @description Reconstruct the cluster variable.
#' Similar to \code{nlme::getGroups}.
#' @name getGroup2-internal
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param data dataset.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param ... [internal] Only used by the generic method.
#'  
#' @return A list containing:
#' \itemize{
#' \item index.cluster: the cluster index for each observation.
#' \item name.cluster: a unique identifier for each cluster.
#' \item n.cluster: the number of clusters.
#' }
#' 
#' @concept extractor
#' @keywords internal
`.getGroups2` <-
    function(object, data, index.Omega, endogenous) UseMethod(".getGroups2")

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
#' .getGroups2(e.lm, data = dW)
#'
#' #### gls model ####
#' e.gls1 <- gls(Y1~X1, data = dW)
#' .getGroups2(e.gls1, data = dW)
#' 
#' e.gls2 <- gls(Y~X1, correlation = corCompSymm(form=~1|id), data = dL)
#' .getGroups2(e.gls2, data = dL)
#'
#' e.gls3 <- gls(Y~X1, correlation = corCompSymm(form=~time|id), data = dL)
#' .getGroups2(e.gls3, data = dL)
#'
#' e.gls4 <- gls(Y~X1, weight = varIdent(form=~1|time2), data = dL)
#' .getGroups2(e.gls4, data = dL)
#'
#' e.gls5 <- gls(Y~X1, weight = varIdent(form=~1|time2),
#'               correlation = corSymm(form=~time|id), data = dL)
#' .getGroups2(e.gls5, data = dL)
#'
#' #### lme model ####
#' e.lme <- lme(Y~X1, random=~1|id, data = dL)
#' .getGroups2(e.lme, data = dL)
#' 
#' #### lvm model ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' .getGroups2(e.lvm, data = dW)

## * .getGroups2.lm
#' @rdname getGroups2-internal
.getGroups2.lm <- function(object, data = NULL, index.Omega = NULL, endogenous){
    if(is.null(data)){
        data <- extractData(object)
    }
    if(is.null(index.Omega)){
        index.Omega <- .getIndexOmega(object, data = data)
    }
    n.obs <- NROW(data)
    out <- list(index.cluster = 1:n.obs,
                name.cluster = 1:n.obs,                
                n.cluster = n.obs,
                index.Omega = as.list(index.Omega)
                )
    return(out)
    
}
    
## * .getGroups2.gls
#' @rdname getGroups2-internal
.getGroups2.gls <- function(object, data = NULL, index.Omega = NULL, endogenous = NULL){
    if(is.null(data)){
        data <- extractData(object)
    }
    if(is.null(index.Omega)){
        index.Omega <- .getIndexOmega(object, data = data)
    }
    if(is.null(endogenous)){
        endogenous <- lava::endogenous(object)
    }
    n.obs <- NROW(data)

    ## ** find clusters
    if(!is.null(object$modelStruct$reStruct)){
        iFormula <- formula(object$modelStruct$reStruct)
        if(length(iFormula)>1){
            stop(".getGroups2 does not handle multiple random effects \n")
        }
        index.cluster <- as.numeric(nlme::getGroups(data,
                                                    form = iFormula[[1]])
                                    )
    }else if(!is.null(object$modelStruct$corStruct)){
        iFormula <- formula(object$modelStruct$corStruct)
        index.cluster <- as.numeric(nlme::getGroups(data,
                                                    form = iFormula)
                                    )
        ## as.numeric(nlme::getGroups(object))
    }else{
        index.cluster <- 1:n.obs
    }
    name.cluster <- unique(index.cluster)
    n.cluster <- length(name.cluster)

    ## ** reorder cluster according to the data ordering
    index.cluster <- as.numeric(factor(index.cluster, levels = name.cluster))

    ## ** export
    return(list(index.cluster = index.cluster,
                name.cluster = name.cluster,
                n.cluster = n.cluster,
                index.Omega = tapply(index.Omega,index.cluster, list)
                ))
}

## * .getGroups2.lme
#' @name getGroups2-internal
.getGroups2.lme <- .getGroups2.gls


## * .getGroups2.lvmfit
#' @rdname getGroups2-internal
.getGroups2.lvmfit <- function(object, data = NULL, index.Omega = NULL, endogenous = NULL){
    if(is.null(data)){
        data <- extractData(object)
    }
    if(is.null(index.Omega)){
        index.Omega <- .getIndexOmega(object, data = data)
    }
    if(is.null(endogenous)){
        endogenous <- lava::endogenous(object)
    }
    n.endogenous <- length(endogenous)

    ## ** find clusters
    n.cluster <- NROW(data)
    name.cluster <- 1:n.cluster
    missing <- any(is.na(index.Omega))
    index.cluster <- unlist(lapply(name.cluster, rep, times = n.endogenous))

    if(missing){
        index.cluster <- index.cluster[!is.na(index.Omega)]
        index.Omega <- tapply(index.Omega[!is.na(index.Omega)], index.cluster, list)
    }else{
        index.Omega <- tapply(index.Omega, index.cluster, list)
    }
    return(list(index.cluster = index.cluster,
                name.cluster = name.cluster,
                n.cluster = n.cluster,
                index.Omega = index.Omega
                ))
    
}










######################################################################
### sCorrect-getGroups2.R ends here
