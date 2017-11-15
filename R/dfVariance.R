### dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: nov 15 2017 (17:46) 
##           By: Brice Ozenne
##     Update #: 137
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - dfVariance
#' @title  Compute the degree of freedom of the variance parameters
#' @description Compute the degree of freedom of the variance parameters
#' @name dfVariance
#' @export
`dfVariance` <-
  function(object, ...) UseMethod("dfVariance")

## * dfVariance.gls
#' @rdname dfVariance
#' @export
dfVariance.gls <- function(object, cluster, ...){

    p <- .coef2(object)
    data <- getData(object)
    vcov.param <- attr(residuals2(object, cluster = cluster, p = p, data = data,
                                  adjust.residuals = FALSE, return.vcov.param = TRUE),
                       "vcov.param")

    N.param <- length(p)
    name.param <- names(p)
    
### ** Define function to compute the standard errors
    calcSigma <- function(iParam){ # x <- p.obj
        M <- attr(residuals2(object, cluster = cluster, p = iParam, data = data,
                             adjust.residuals = FALSE, return.vcov.param = TRUE),
                  "vcov.param")
        return(setNames(diag(M), name.param))         
    }

### ** Compute the gradient of the function computing the standard errors
    browser()
    dSigma.dtheta <- numDeriv::jacobian(func = calcSigma, x = p, method = "Richardson")
    
### ** Compute degrees of freedom

    ## diag(vcov.param) - calcSigma(p)
    numerator <- 2*diag(vcov.param)^2
    denom <- rowSums(dSigma.dtheta %*% vcov.param * dSigma.dtheta)
    df <- numerator/denom
    alpha <- df/diag(vcov.param)
    return(df)
 
}

## * dfVariance.lme
#' @rdname dfVariance
#' @export
dfVariance.lme <- dfVariance.gls

## * dfVariance.lvmfit
#' @rdname dfVariance
#' @export
dfVariance.lvmfit <- function(object, p = NULL, data = NULL, vcov.param = NULL,
                              robust = FALSE, adjust.residuals = FALSE, power = 0.5, as.clubSandwich = TRUE, ...){

    if(is.null(p)){p <- pars(object)}
    if(is.null(data)){data <- model.frame(object)}
    if(is.null(vcov.param)){
        if(robust==FALSE){
            vcov.param <- vcov(object)
        }else{
            vcov.param <- crossprod(iid2(object, p = p, data = data,
                                         adjust.residuals = adjust.residuals,
                                         power = power,
                                         as.clubSandwich = as.clubSandwich))
        }
    }

    n.param <- length(p)
    name.param <- names(p)
    
### ** Define function to compute the standard errors
    if(robust == FALSE){
        if(adjust.residuals==FALSE){
            calcSigma <- function(iParam){ # x <- p.obj
                M <- solve(information(object, p = iParam))
                return(setNames(diag(M), name.param))         
            }
        }else{
            stop("adjust.residuals=TRUE not implemented for standard errors derived obtained from the information matrix \n")
        }
    }else{
        calcSigma <- function(iParam){ # x <- p.obj
            M <- crossprod(iid2(object, p = iParam, data = data,
                                adjust.residuals = adjust.residuals,
                                power = power,
                                as.clubSandwich = as.clubSandwich))
            return(setNames(diag(M), name.param))         
        }
    }    

### ** Compute the gradient of the function computing the standard errors
    dSigma.dtheta <- numDeriv::jacobian(func = calcSigma, x = p, method = "Richardson")
    
### ** Compute degrees of freedom

    ## diag(vcov.param) - calcSigma(p)
    numerator <- 2*diag(vcov.param)^2
    denom <- rowSums(dSigma.dtheta %*% vcov.param * dSigma.dtheta)
    df <- numerator/denom
    alpha <- df/diag(vcov.param)
    return(df)
}



##----------------------------------------------------------------------
### dfVariance.R ends here
