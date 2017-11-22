### dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: nov 17 2017 (10:33) 
##           By: Brice Ozenne
##     Update #: 219
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
#'
#' @param object a lvm object
#' @param cluster the grouping variable relative to which the observations are i.i.d.
#' 
#' @export
`dfVariance` <-
  function(object, ...) UseMethod("dfVariance")

## * dfVariance.gls
#' @rdname dfVariance
#' @export
dfVariance.gls <- function(object, cluster,
                           robust = FALSE, adjust.residuals = FALSE, power = 0.5, as.clubSandwich = TRUE, ...){

    p <- .coef2(object)
    data <- getData(object)
    N.param <- length(p)
    name.param <- names(p)
    
### ** Define function to compute the standard errors
    if(robust){
        stop("not implemented \n")
    }else{
        calcSigma <- function(iParam){ # x <- p.obj
            pp <- p
            pp[names(iParam)] <- iParam
            return(attr(residuals2(object, cluster = cluster, p = pp, data = data,
                                 adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich,
                                 return.vcov.param = TRUE),
                      "vcov.param"))
        }
    }
    
    calcDiagSigma <- function(iParam){
        return(setNames(diag(calcSigma(iParam)), name.param))         
    }

### ** Compute the gradient of the function computing the standard errors
    keep.param <- setdiff(name.param,attr(.coef2(object),"mean.coef"))
    jac.param <- p[keep.param]
    dSigma.dtheta <- numDeriv::jacobian(func = calcDiagSigma, x = jac.param, method = "Richardson")
    
### ** Compute degrees of freedom

    ## diag(vcov.param) - calcSigma(p)
    vcov.param <- calcSigma(p)
    numerator <- 2*diag(vcov.param)^2
    denom <- rowSums(dSigma.dtheta %*% vcov.param[keep.param,keep.param,drop=FALSE] * dSigma.dtheta)
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
dfVariance.lvmfit <- function(object, fix.mean = FALSE, 
                              robust = FALSE, adjust.residuals = FALSE, power = 0.5, as.clubSandwich = TRUE, ...){

    p <- pars(object)
    data <- model.frame(object)
 
    n.param <- length(p)
    name.param <- names(p)
    
### ** Define function to compute the information matrix
    if(robust == FALSE){
        if(adjust.residuals==FALSE){
            calcI <- function(iParam){ # x <- p.obj
                pp <- p
                pp[names(iParam)] <- iParam
                return(information(object, p = pp))
            }
        }else{
            calcVcov <- function(iParam){ # x <- p.obj
                pp <- p
                pp[names(iParam)] <- iParam                
                return(attr(residuals2(object, p = pp, data = data,
                                              adjust.residuals = adjust.residuals,
                                              power = power,
                                              as.clubSandwich = as.clubSandwich,
                                              return.vcov.param = TRUE), "vcov.param"))
            }
        }
    }else{
        calcVcov <- function(iParam){ # x <- p.obj
            pp <- p
            pp[names(iParam)] <- iParam
            return(crossprod(iid2(object, p = pp, data = data,
                                  adjust.residuals = adjust.residuals,
                                  power = power,
                                  as.clubSandwich = as.clubSandwich)))
        }
        
    }
    
### ** Compute variance-covariance matrix
    if(robust == FALSE && adjust.residuals==FALSE){
        vcov.param <- chol2inv(chol(calcI(p)))
    }else{
        vcov.param <- calcVcov(p)
    }
    dimnames(vcov.param) <- list(name.param, name.param)
    
### ** Compute the gradient of the function computing the standard errors
    if(robust == FALSE && adjust.residuals==FALSE){
        calcDiagVcov <- function(iParam){
            dVcov.dtheta <-  - vcov.param %*% calcI(iParam) %*% vcov.param
            return(as.double(diag(dVcov.dtheta)))         
        }
    }else{
        calcDiagVcov <- function(iParam){
            return(as.double(diag(calcVcov(iParam))))         
        }
    }
    
    keep.type <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
    tableType <- coefType(object, as.lava=FALSE)        
    keep.param <- tableType[!is.na(lava) & detail%in%keep.type,originalLink]

    jac.param <- p[keep.param]
    dVcov.dtheta <- numDeriv::jacobian(calcDiagVcov, x = jac.param, method = "Richardson")

### ** Compute degrees of freedom

    ## diag(vcov.param) - calcVcov(p)
    numerator <- 2*diag(vcov.param)^2
    denom <- rowSums(dVcov.dtheta %*% vcov.param[keep.param,keep.param,drop=FALSE] * dVcov.dtheta)
    df <- numerator/denom
    alpha <- df/diag(vcov.param)
    return(df)
}



##----------------------------------------------------------------------
### dfVariance.R ends here
