### sCorrect-getGroups2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:58) 
## Version: 
## Last-Updated: Jan  3 2022 (09:51) 
##           By: Brice Ozenne
##     Update #: 154
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
#' #### latent variable model ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' .getGroups2(e.lvm, data = dW)

## * .getGroups2.lvm
#' @rdname getGroups2-internal
.getGroups2.lvm <- function(object, data = NULL, index.Omega = NULL, endogenous = NULL){
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
    Uindex.Omega <- unique(index.Omega)
    return(list(index.cluster = index.cluster,
                name.cluster = name.cluster,
                n.cluster = n.cluster,
                index.Omega = index.Omega,
                index2endogenous = setNames(as.list(Uindex.Omega),Uindex.Omega)
                ))
    
}

## * .getGroups2.lvmfit
#' @rdname getGroups2-internal
.getGroups2.lvmfit <- .getGroups2.lvm








######################################################################
### sCorrect-getGroups2.R ends here
