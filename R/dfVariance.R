### dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: nov 30 2017 (15:28) 
##           By: Brice Ozenne
##     Update #: 274
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

## * dfVariance.lm
#' @rdname dfVariance
#' @export
dfVariance.lm <- function(object, adjust.residuals = FALSE, ...){
    name.coef <- names(coef(object))
    n.coef <- length(name.coef)
    df <- setNames(rep(NA,n.coef), name.coef)

    n <- NROW(object$model)
    p <- object$rank

    if(adjust.residuals){
        df[] <- n*(n/(n+p))^2
    }else{
        df[] <- n
    }
    return(df)
}
     

## * dfVariance.gls
#' @rdname dfVariance
#' @export
dfVariance.gls <- function(object, cluster,
                           robust = FALSE, adjust.residuals = FALSE, ...){

    p <- .coef2(object)
    data <- getData(object)
    N.param <- length(p)
    name.param <- names(p)

    power <- 0.5
    as.clubSandwich <- 2
    
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
dfVariance.lvmfit <- function(object, 
                              robust = FALSE, adjust.residuals = FALSE, ...){

    p <- pars(object)
    data <- model.frame(object)
 
    n.param <- length(p)
    name.param <- names(p)

    power <- 0.5
    as.clubSandwich <- 2
    
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

     vcov.param/vcov(e)
    browser()
    vcov.param[1:4,1:4] - coef(e)["Y~~Y"] * solve(t(X) %*% X)
    dVcov.dtheta[1:4] - diag(vcov.param)[1:4] / coef(e)["Y~~Y"]

    diag(vcov.param)[1:4]^3 / coef(e)["Y~~Y"]^2
    vcov.param[keep.param,keep.param,drop=FALSE] -  (2/n*coef(e)["Y~~Y"]^(2))
### ** Compute degrees of freedom

    ## diag(vcov.param) - calcVcov(p)
    numerator <- 2*diag(vcov.param)^2
    denom <- rowSums(dVcov.dtheta %*% vcov.param[keep.param,keep.param,drop=FALSE] * dVcov.dtheta)
    df <- numerator/denom
    alpha <- df/diag(vcov.param)
    return(df)
}


## * .dinformation2
.dinformation2 <- function(dmu.dtheta, d2mu.d2theta,
                           dOmega.dtheta, d2Omega.d2theta,
                           Omega, ls.indexOmega, bias.Omega,
                           n.param, name.param, n.cluster){

### ** prepare
    clusterSpecific <- !is.null(ls.indexOmega)
    correction <- !is.null(bias.Omega)
    if(clusterSpecific==FALSE){

        if(correction){
            Tbias.vcov.param <- Reduce("+",bias.Omega)/n.cluster
            Omega <- Omega + Tbias.vcov.param
        }
        
        iOmega <- chol2inv(chol(Omega))
        
    }

### ** compute the derivative of the information matrix for each parameters
    dInfo <-  matrix(0, nrow = n.param, ncol = n.param, dimnames = list(name.param,name.param))
    
    for(iP1 in 1:n.param){ # iP <- 1
        for(iP3 in 1:n.param){ # iP <- 1
                
                iName1 <- name.param[iP1]
                iName3 <- name.param[iP3]

                test.Omega <- !is.null(dOmega.dtheta[[iName1]]) && !is.null(dOmega.dtheta[[iName3]])
                test.mu <- !is.null(dmu.dtheta[[iName1]])
                
                ## *** Individual specific Omega (e.g. presence of missing values)
                if(clusterSpecific){
                    for(iC in 1:n.cluster){
                    
                    Omega.tempo <- Omega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                    if(correction){
                        Omega.tempo <- Omega.tempo + bias.Omega[[iC]]
                    }
                    iOmega.tempo <- chol2inv(chol(Omega.tempo))

                    if(test.Omega){
                        iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        iOmega.dOmega.3 <- iOmega.tempo %*% dOmega.dtheta[[iName3]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        
                        dInfo[iP1,iP3] <- dInfo[iP1,iP3] + tr(iOmega.dOmega.3 %*% iOmega.dOmega.1 %*% iOmega.dOmega.1)
                        
                        if(!is.null(d2Omega.d2theta[[iName1]][[iName3]])){
                            iOmega.d2Omega.13 <- iOmega.tempo %*% d2Omega.d2theta[[iName1]][[iName3]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                            dInfo[iP1,iP3] <- dInfo[iP1,iP3] + tr(iOmega.d2Omega.13 %*% iOmega.dOmega.1)                
                        }else if(!is.null(d2Omega.d2theta[[iName3]][[iName1]])){
                            iOmega.d2Omega.13 <- iOmega.tempo %*% d2Omega.d2theta[[iName3]][[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                            dInfo[iP1,iP3] <- dInfo[iP1,iP3] + tr(iOmega.d2Omega.13 %*% iOmega.dOmega.1)                
                        }                        
                    }                    
                    
                    if(test.mu){
                        dmu.1 <- dmu.dtheta[[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]

                        if(!is.null(d2mu.d2theta[[name1]][[name3]])){
                            dInfo[iP1,iP3] <- dInfo[iP1,iP3] + 2 * sum(d2mu.d2theta[[name1]][[name3]] %*% iOmega.tempo * dmu.1)
                        }else if(!is.null(d2mu.d2theta[[name3]][[name1]])){
                            dInfo[iP1,iP3] <- dInfo[iP1,iP3] + 2 * sum(d2mu.d2theta[[name3]][[name1]] %*% iOmega.tempo * dmu.1)
                        }
                        if(!is.null(dOmega.dtheta[[name3]])){
                            dOmega.3 <- dOmega.dtheta[[iName3]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                            dInfo[iP1,iP3] <- dInfo[iP1,iP3] + sum(dmu.tempo1 %*% iOmega.tempo %*% dOmega.3 %*% iOmega.tempo * dmu.tempo1)
                        }
                        
                    }
                    
                }
            }
            
            ## *** Same for all individuals
            if(clusterSpecific == FALSE){
                if(test.Omega){
                    iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]]
                    iOmega.dOmega.3 <- iOmega.tempo %*% dOmega.dtheta[[iName3]]

                    dInfo[iP1,iP3] <- dInfo[iP1,iP3] + n.cluster * tr(iOmega.dOmega.3 %*% iOmega.dOmega.1 %*% iOmega.dOmega.1)
                    
                    if(!is.null(d2Omega.d2theta[[iName1]][[iName3]])){
                        dInfo[iP1,iP3] <- dInfo[iP1,iP3] + n.cluster * tr(iOmega.tempo %*% d2Omega.d2theta[[iName1]][[iName3]] %*% iOmega.dOmega.1)                
                    }else if(!is.null(d2Omega.d2theta[[iName3]][[iName1]])){
                        dInfo[iP1,iP3] <- dInfo[iP1,iP3] + n.cluster * tr(iOmega.tempo %*% d2Omega.d2theta[[iName3]][[iName1]] %*% iOmega.dOmega.1)                
                    }
                    
                }
                                
                if(test.mu){

                    if(!is.null(d2mu.d2theta[[name1]][[name3]])){
                        dInfo[iP1,iP3] <- dInfo[iP1,iP3] + 2 * sum(d2mu.d2theta[[name1]][[name3]] %*% iOmega.tempo * dmu.dtheta[[iName1]])
                    }else if(!is.null(d2mu.d2theta[[name3]][[name1]])){
                        dInfo[iP1,iP3] <- dInfo[iP1,iP3] + 2 * sum(d2mu.d2theta[[name3]][[name1]] %*% iOmega.tempo * dmu.dtheta[[iName1]])
                    }

                    if(!is.null(dOmega.dtheta[[name3]])){
                        dInfo[iP1,iP3] <- dInfo[iP1,iP3] + sum(dmu.dtheta[[iName1]] %*% iOmega.tempo %*% dOmega.dtheta[[name3]] %*% iOmega.tempo * dmu.dtheta[[iName1]])
                    }

                }

            }
        }
    }

### ** export
    return(dInfo)
}


##----------------------------------------------------------------------
### dfVariance.R ends here
