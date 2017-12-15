### dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: dec 15 2017 (16:55) 
##           By: Brice Ozenne
##     Update #: 468
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
dfVariance.lm <- function(object, adjust.residuals = TRUE, ...){
    object.coef <- coef(object)
    name.coef <- names(object.coef)
    n.coef <- length(name.coef)
    df <- setNames(rep(NA,n.coef+1), c(name.coef,"sigma"))

    n <- NROW(object$model)
    p <- object$rank

    if(adjust.residuals==FALSE){
        df[name.coef] <- n
        df["sigma"] <- n/4
    }else{
        df[name.coef] <- n^2/(n+p)
        df["sigma"] <- n^2/(4*(n+p))
    }
    
    return(df)
}
     

## * dfVariance.gls
#' @rdname dfVariance
#' @export
dfVariance.gls <- function(object, cluster, vcov.param = NULL,
                           adjust.residuals = FALSE, numericDerivative = FALSE, ...){

    p <- .coef2(object)
    data <- getData(object)
    N.param <- length(p)
    name.param <- names(p)

    power <- 0.5
    as.clubSandwich <- 1

    keep.param <- setdiff(name.param, attr(.coef2(object),"mean.coef"))

    ### ** Compute the covariance matrix
    calcSigma <- function(iParam){ # x <- p.obj
        pp <- p
        pp[names(iParam)] <- iParam
        return(attr(residuals2(object, cluster = cluster, p = pp, data = data,
                               adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich,
                               return.vcov.param = TRUE, second.order = FALSE),
                    "vcov.param"))
    }
    if(is.null(vcov.param)){
        vcov.param <- calcSigma(p)
    }
    
    ### ** Compute the gradient of the standard errors
    if(numericDerivative){
    
        calcDiagSigma <- function(iParam){
            return(setNames(diag(calcSigma(iParam)), name.param))         
        }

        jac.param <- p[keep.param]
        dSigma.dtheta <- numDeriv::jacobian(func = calcDiagSigma, x = jac.param, method = "Richardson")
        
    }else{

        dots <- list(...)
        if("prepareScore" %in% names(dots)){
            prepareScore <- dots$prepareScore
        }else{
            prepareScore  <- attr(residuals2(object, cluster = cluster, p = p, data = data,
                                        adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich,
                                        return.prepareScore2 = TRUE, second.order = TRUE), "prepareScore2")
        }

        dInfo.dtheta <- .dinformation2(dmu.dtheta = prepareScore$dmu.dtheta,
                                       d2mu.d2theta = NULL,
                                       dOmega.dtheta = prepareScore$dOmega.dtheta,
                                       d2Omega.d2theta = prepareScore$d2Omega.d2theta,
                                       Omega = prepareScore$Omega,
                                       ls.indexOmega = prepareScore$ls.indexOmega,
                                       hat = prepareScore$hat,
                                       n.param  = prepareScore$n.param,
                                       name.param  = prepareScore$name.param,
                                       name.deriv = keep.param,
                                       n.cluster = prepareScore$n.cluster)

        if(dim(dInfo.dtheta)[3]==1){
            dSigma.dtheta <- -diag(vcov.param %*% dInfo.dtheta[,,1] %*% vcov.param)
        }else{
            dSigma.dtheta <- apply(dInfo.dtheta, 3, function(x){ - diag(vcov.param %*% x %*% vcov.param)})
        }
    }
    print(dSigma.dtheta)
    
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
dfVariance.lvmfit <- function(object, adjust.residuals = FALSE, ...){

    p <- pars(object)
    data <- model.frame(object)
 
    n.param <- length(p)
    name.param <- names(p)

    power <- 0.5
    as.clubSandwich <- 1
    
    ### ** Define function to compute the information matrix
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
    
### ** Compute variance-covariance matrix
    if(adjust.residuals==FALSE){
        vcov.param <- chol2inv(chol(calcI(p)))
    }else{
        vcov.param <- calcVcov(p)
    }    
    dimnames(vcov.param) <- list(name.param, name.param)
    
### ** Compute the gradient of the function computing the standard errors
    if(adjust.residuals==FALSE){
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
    dimnames(dVcov.dtheta) <- list(names(p), names(jac.param))
    
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
                           Omega, ls.indexOmega, hat,
                           n.param, name.param, name.deriv,
                           n.cluster){

### ** prepare
    index.deriv <- match(name.deriv, name.param)
    clusterSpecific <- !is.null(ls.indexOmega)
    iOmega <- chol2inv(chol(Omega))        

    if(!clusterSpecific){## small sample correction                    
        df.mean <- Reduce("+",hat)
        iN.cluster <- as.double(n.cluster - diag(df.mean))
    }
    
    ### ** compute the derivative of the information matrix for each parameters
    dInfo <-  array(0,
                    dim = c(n.param, n.param, length(name.deriv)),
                    dimnames = list(name.param, name.param, name.deriv))
    
    for(iDeriv in index.deriv){ # iDeriv <- 4
        for(iP1 in 1:n.param){ # iP1 <- 1
            for(iP2 in iP1:n.param){ # iP2 <- 1
                
                iNameD <- name.param[iDeriv]
                iName1 <- name.param[iP1]
                iName2 <- name.param[iP2]
                
                test.Omega1 <- !is.null(dOmega.dtheta[[iNameD]]) && !is.null(dOmega.dtheta[[iName1]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega2a <- !is.null(d2Omega.d2theta[[iNameD]][[iName1]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega2b <- !is.null(d2Omega.d2theta[[iName1]][[iNameD]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega3a <- !is.null(d2Omega.d2theta[[iNameD]][[iName2]]) && !is.null(dOmega.dtheta[[iName1]])
                test.Omega3b <- !is.null(d2Omega.d2theta[[iName2]][[iNameD]]) && !is.null(dOmega.dtheta[[iName1]])
                
                test.mu1a <- !is.null(d2mu.d2theta[[iNameD]][[iName1]]) && !is.null(dmu.dtheta[[iName2]])
                test.mu1b <- !is.null(d2mu.d2theta[[iName1]][[iNameD]]) && !is.null(dmu.dtheta[[iName2]])
                test.mu2a <- !is.null(d2mu.d2theta[[iNameD]][[iName2]]) && !is.null(dmu.dtheta[[iName1]])
                test.mu2b <- !is.null(d2mu.d2theta[[iName2]][[iNameD]]) && !is.null(dmu.dtheta[[iName1]])
                test.mu3 <- !is.null(dOmega.dtheta[[iNameD]]) && !is.null(dmu.dtheta[[iName1]]) && !is.null(dmu.dtheta[[iName2]])

                ## *** Individual specific Omega (e.g. presence of missing values)
                if(clusterSpecific){
                    
                    for(iC in 1:n.cluster){

                        ## prepare
                        Omega.tempo <- Omega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        iOmega.tempo <- iOmega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        if(!is.null(dmu.dtheta[[iName1]])){
                            dmu.1 <- dmu.dtheta[[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dmu.dtheta[[iName2]])){
                            dmu.2 <- dmu.dtheta[[iName2]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.mu1a){
                            d2mu.D1 <- d2mu.d2theta[[iNameD]][[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.mu1b){
                            d2mu.D1 <- d2mu.d2theta[[iName1]][[iNameD]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.mu1a){
                            d2mu.D <- d2mu.d2theta[[iNameD]][[iName2]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.mu1b){
                            d2mu.D <- d2mu.d2theta[[iName2]][[iNameD]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iNameD]])){
                            iOmega.dOmega.D <- iOmega.tempo %*% dOmega.dtheta[[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iName1]])){
                            iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iName2]])){
                            iOmega.dOmega.2 <- iOmega.tempo %*% dOmega.dtheta[[iName2]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.Omega2a){
                            d2Omega.D1 <- d2Omega.d2theta[[iNameD]][[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.Omega2b){
                            d2Omega.D1 <- d2Omega.d2theta[[iName1]][[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.Omega3a){
                            d2Omega.D2 <- d2Omega.d2theta[[iNameD]][[iName2]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }else{
                            d2Omega.D2 <- d2Omega.d2theta[[iName2]][[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }

                        ## small sample correction  
                        iW.cluster <- 1 -  diag(hat[[iC]])
                        
                        ## compute
                        if(test.Omega1){
                            iDiag1 <- diag(iOmega.dOmega.D %*% iOmega.dOmega.1 %*% iOmega.dOmega.2)
                            iDiag2 <- diag(iOmega.dOmega.1 %*% iOmega.dOmega.D %*% iOmega.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iW.cluster + iDiag2 * iW.cluster)
                        }
                        
                        if(test.Omega2a || test.Omega2b){                            
                            iDiag <- diag(iOmega %*% d2Omega.D1 %*% iOmega.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iW.cluster)
                        }

                        if(test.Omega3a || test.Omega3b){
                            iDiag <- diag(iOmega.dOmega.1 %*% iOmega %*% d2Omega.D2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iW.cluster)
                        }

                        if(test.mu1a || test.mu1b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + d2mu.D1 %*% iOmega %*% t(dmu.2)
                        }

                        if(test.mu2a || test.mu2b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + dmu.1 %*% iOmega %*% t(d2mu.D2)
                        }

                        if(test.mu3){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - dmu.1 %*% iOmega.dOmega.D %*% iOmega.tempo %*% t(dmu.2)
                        }
                    
                }
            }
            
            ## *** Same for all individuals
                if(clusterSpecific == FALSE){

                    ## prepare
                    if(!is.null(dmu.dtheta[[iName1]])){
                        dmu.1 <- dmu.dtheta[[iName1]]
                    }
                    if(!is.null(dmu.dtheta[[iName2]])){
                        dmu.2 <- dmu.dtheta[[iName2]]
                    }
                    if(test.mu1a){
                        d2mu.D1 <- d2mu.d2theta[[iNameD]][[iName1]]
                    }else if(test.mu1b){
                        d2mu.D1 <- d2mu.d2theta[[iName1]][[iNameD]]
                    }
                    if(test.mu1a){
                        d2mu.D2 <- d2mu.d2theta[[iNameD]][[iName2]]
                    }else if(test.mu1b){
                        d2mu.D2 <- d2mu.d2theta[[iName2]][[iNameD]]
                    }
                    if(!is.null(dOmega.dtheta[[iNameD]])){
                        iOmega.dOmega.D <- iOmega.tempo %*% dOmega.dtheta[[iNameD]]
                    }
                    if(!is.null(dOmega.dtheta[[iName1]])){
                        iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]]
                    }
                    if(!is.null(dOmega.dtheta[[iName2]])){
                        iOmega.dOmega.2 <- iOmega.tempo %*% dOmega.dtheta[[iName2]]
                    }
                    if(test.Omega2a){
                        d2Omega.D1 <- d2Omega.d2theta[[iNameD]][[iName1]]
                    }else if(test.Omega2b){
                        d2Omega.D1 <- d2Omega.d2theta[[iName1]][[iNameD]]
                    }
                    if(test.Omega3a){
                        d2Omega.D2 <- d2Omega.d2theta[[iNameD]][[iName2]]
                    }else{
                        d2Omega.D2 <- d2Omega.d2theta[[iName2]][[iNameD]]
                    }

                    ## compute
                    if(test.Omega1){
                        iDiag1 <- diag(iOmega.dOmega.D %*% iOmega.dOmega.1 %*% iOmega.dOmega.2)
                        iDiag2 <- diag(iOmega.dOmega.1 %*% iOmega.dOmega.D %*% iOmega.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iN.cluster + iDiag2 * iN.cluster)
                    }

                    if(test.Omega2a || test.Omega2b){
                        iDiag <- diag(iOmega %*% d2Omega.D1 %*% iOmega.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.cluster)
                    }

                    if(test.Omega3a || test.Omega3b){
                        iDiag <- diag(iOmega.dOmega.1 %*% iOmega %*% d2Omega.D2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.cluster)
                    }

                    if(test.mu1a || test.mu1b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + d2mu.D1 %*% iOmega %*% t(dmu.2)
                    }

                    if(test.mu2a || test.mu2b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + dmu.1 %*% iOmega %*% t(d2mu.D2)
                    }

                    if(test.mu3){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - dmu.1 %*% iOmega.dOmega.D %*% iOmega.tempo %*% t(dmu.2)
                    }
            }
            
            }
        }
        dInfo[,,iNameD] <- dInfo[,,iNameD] + t(dInfo[,,iNameD]) - diag(diag(dInfo[,,iNameD]))
    }

    ### ** export
    return(dInfo)
}


##----------------------------------------------------------------------
### dfVariance.R ends here
