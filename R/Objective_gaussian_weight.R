### Objective_gaussian_weight.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 17 2020 (16:29) 
## Version: 
## Last-Updated: feb 21 2020 (09:48) 
##           By: Brice Ozenne
##     Update #: 89
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Estimate LVM With Weights
##' @description Estimate LVM with weights.
##' @name gaussian_weight
##' 
##' @examples
##' #### linear regression with weights ####
##'
##' ## data
##' df <- data.frame(Y = c(1,2,2,1,2),
##'                  X = c(1,1,2,2,2),
##'                  missing = c(0,0,0,0,1),
##'                  weights = c(1,1,2,1,NA))
##'
##' ## using lm
##' e.lm.GS <- lm(Y~X, data = df)
##' e.lm.test <- lm(Y~X, data = df[df$missing==0,], weights = df[df$missing==0,"weights"])
##' 
##' ## using lvm
##' m <- lvm(Y~X)
##' e.GS <- estimate(m, df)
##' e.test <- estimate(m, data = df[df$missing==0,],
##'                    weights = df[df$missing==0,"weights"],
##'                    estimator = "gaussian_weight")
##' 

## * gaussian_weight.estimate.hook
##' @name gaussian_weight
##' @export
gaussian_weight.estimate.hook <- function(x, data, estimator, ...){
    dots <- list(...)
    if(identical(estimator,"gaussian_weight")){
        xe <- suppressWarnings(estimate(x, data = data, control = list(iter.max = 0))) ## initialize coefficients
        x$sCorrect <- conditionalMoment(xe, data = data, param = initCoef,
                                        initialize = TRUE, first.order = TRUE, second.order = FALSE, usefit = FALSE)
    }
    return( c(list(x=x, data=data, estimator = estimator),dots) )
}

##' @name gaussian_weight
##' @export
gaussian_weight_method.lvm <- "nlminb2"

## * gaussian_weight_logLik.lvm
##' @name gaussian_weight
##' @export
`gaussian_weight_logLik.lvm` <- function(object,type="cond",p,data,weights,...) {
    ## ** compute mu and Omega
    if(type!="cond"){
        stop("Not implemented for other types than \"cond\"\n ")
    }
    cM <- conditionalMoment(object,
                            initialize = FALSE, first.order = FALSE, second.order = FALSE,
                            usefit = TRUE, residuals = TRUE, leverage = FALSE,
                            param = p,
                            data = data)
    
    ## ** prepare
    name.pattern <- cM$missing$name.pattern
    missing.pattern <- cM$missing$pattern
    unique.pattern <- cM$missing$unique.pattern
    n.pattern <- length(name.pattern)
    
    OmegaM1 <- cM$moment$OmegaM1
    residuals <- cM$residuals
    logLik <- 0

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)
        iResiduals <- residuals[iIndex,iY,drop=FALSE]
        iM <- length(iY)
        if(is.null(weights)){
            logLik <- logLik - (cM$cluster$n.cluster/2) * (iM * log(2*pi) - log(det(iOmegaM1))) - sum((iResiduals %*% iOmegaM1) * iResiduals)/2
        }else{
            logLik <- logLik - (sum(weights)/2) * (iM * log(2*pi) - log(det(iOmegaM1))) - sum(weights[,1]/2 * rowSums((iResiduals %*% iOmegaM1) * iResiduals))
        }
    }
    return(logLik)
}

##' @name gaussian_weight
##' @export
`gaussian_weight_objective.lvm` <- function(x, ...) {
    logLik <- gaussian_weight_logLik.lvm(object = x,...)
    return(-logLik)
}

## * gaussian_weight_score.lvm
##' @name gaussian_weight
##' @export
gaussian_weight_score.lvm <- function(x, data, p, S, n, mu=NULL, weights=NULL, debug=FALSE, reindex=FALSE, mean=TRUE, constrain=TRUE, indiv=FALSE,...) {

    ## if(constrain){
    ##     stop("gaussian_weight_score.lvm does not handle constrain")
    ## }
    ## if(reindex){
    ##     stop("gaussian_weight_score.lvm does not handle reindex")
    ## }
    ## if(!mean){
    ##     stop("gaussian_weight_score.lvm only handles mean")
    ## }
    
    ## ** compute moments
    cM <- conditionalMoment(x,
                            initialize = FALSE, first.order = TRUE, second.order = FALSE,
                            usefit = TRUE, residuals = TRUE, leverage = FALSE,
                            param = p,
                            data = data)

    ## ** compute score
    score <- .score2(dmu = cM$dmoment$dmu,
                     dOmega = cM$dmoment$dOmega,                    
                     epsilon = cM$residuals,
                     OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                     missing.pattern = cM$missing$pattern,
                     unique.pattern = cM$missing$unique.pattern,
                     name.pattern = cM$missing$name.pattern,
                     name.param = cM$skeleton$Uparam,
                     name.meanparam = cM$skeleton$Uparam.mean,
                     name.varparam = cM$skeleton$Uparam.var,
                     weights = weights[,1],
                     n.cluster = cM$cluster$n.cluster)

    ## ** export
    if(indiv){
        return(score)
    }else{
        return(colSums(score))
    }
}

##' @name gaussian_weight
##' @export
gaussian_weight_gradient.lvm <-  function(...) {
    return(-gaussian_weight_score.lvm(...))
}

## * gaussian_weight_hessian.lvm
##' @name gaussian_weight
##' @export
`gaussian_weight_hessian.lvm` <- function(x,p,n, weights=NULL,...) {

    ## ** compute moments
    cM <- conditionalMoment(x,
                            initialize = FALSE, first.order = TRUE, second.order = FALSE,
                            usefit = TRUE, residuals = FALSE, leverage = TRUE,
                            param = p,
                            data = data)
    ## ** compute hessian
    info <- .information2(dmu = cM$dmoment$dmu,
                         dOmega = cM$dmoment$dOmega,
                         OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                         missing.pattern = cM$missing$pattern,
                         unique.pattern = cM$missing$unique.pattern,
                         name.pattern = cM$missing$name.pattern,
                         grid.mean = cM$skeleton$grid.dmoment$mean, 
                         grid.var = cM$skeleton$grid.dmoment$var, 
                         name.param = cM$skeleton$Uparam,
                         leverage = cM$leverage,
                         weights = weights[,1],
                         n.cluster = cM$cluster$n.cluster)

    ## ** export
    return(info)
}

######################################################################
### Objective_gaussian_weight.R ends here
