### dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: nov 16 2017 (16:51) 
##           By: Brice Ozenne
##     Update #: 159
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
#' @param fast should the degree of freedom be computed keeping the mean parameters fixed?
#' 
#' @export
`dfVariance` <-
  function(object, ...) UseMethod("dfVariance")

## * dfVariance.gls
#' @rdname dfVariance
#' @export
dfVariance.gls <- function(object, cluster, fast = TRUE, ...){

    p <- .coef2(object)
    data <- getData(object)
    vcov.param <- attr(residuals2(object, cluster = cluster, p = p, data = data,
                                  adjust.residuals = FALSE, return.vcov.param = TRUE),
                       "vcov.param")

    N.param <- length(p)
    name.param <- names(p)
    
### ** Define function to compute the standard errors
    calcSigma <- function(iParam){ # x <- p.obj
        pp <- p
        pp[names(iParam)] <- iParam
        M <- attr(residuals2(object, cluster = cluster, p = pp, data = data,
                             adjust.residuals = FALSE, return.vcov.param = TRUE),
                  "vcov.param")
        return(setNames(diag(M), name.param))         
    }

### ** Compute the gradient of the function computing the standard errors
    if(fast){
        keep.param <- setdiff(name.param,attr(.coef2(object),"mean.coef"))
    }else{
        keep.param <- name.param
    }
    jac.param <- p[keep.param]
    dSigma.dtheta <- numDeriv::jacobian(func = calcSigma, x = jac.param, method = "Richardson")
    
### ** Compute degrees of freedom

    ## diag(vcov.param) - calcSigma(p)
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
dfVariance.lvmfit <- function(object, fast = TRUE, vcov.param = NULL,
                              robust = FALSE, adjust.residuals = FALSE, power = 0.5, as.clubSandwich = TRUE, ...){

    p <- pars(object)
    data <- model.frame(object)
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
                pp <- p
                pp[names(iParam)] <- iParam
                M <- solve(information(object, p = pp))
                return(setNames(diag(M), name.param))         
            }
        }else{
            stop("adjust.residuals=TRUE not implemented for standard errors derived obtained from the information matrix \n")
        }
    }else{
        calcSigma <- function(iParam){ # x <- p.obj
            pp <- p
            pp[names(iParam)] <- iParam
            M <- crossprod(iid2(object, p = pp, data = data,
                                adjust.residuals = adjust.residuals,
                                power = power,
                                as.clubSandwich = as.clubSandwich))
            return(setNames(diag(M), name.param))         
        }
    }    

### ** Compute the gradient of the function computing the standard errors
    if(fast){
        keep.type <- c("Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
        tableType <- coefType(object, as.lava=FALSE)        
        keep.param <- tableType[!is.na(lava) & detail%in%keep.type,originalLink]
    }else{
        keep.param <- name.param
    }
    jac.param <- p[keep.param]
    dSigma.dtheta <- numDeriv::jacobian(func = calcSigma, x = jac.param, method = "Richardson")
    
### ** Compute degrees of freedom

    ## diag(vcov.param) - calcSigma(p)
    numerator <- 2*diag(vcov.param)^2
    denom <- rowSums(dSigma.dtheta %*% vcov.param[keep.param,keep.param,drop=FALSE] * dSigma.dtheta)
    df <- numerator/denom
    alpha <- df/diag(vcov.param)
    return(df)
}



##----------------------------------------------------------------------
### dfVariance.R ends here
