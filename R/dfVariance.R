### dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: jan  3 2018 (18:40) 
##           By: Brice Ozenne
##     Update #: 512
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
dfVariance.lvmfit <- function(object, C = NULL, ...){

    p <- pars(object)
    q <- NROW(C)
    n.param <- length(p)
    name.param <- names(p)

    if(is.null(C)){
        C <- diag(1, nrow = n.param, ncol = n.param)
        dimnames(C) <- list(name.param, name.param)
    }else{
        if(NCOL(C) != n.param){
            stop("Argument \'C\' should be a matrix with ",n.param," columns \n")
        }
        if(any(abs(svd(C)$d)<1e-10)){
            stop("Argument \'C\' is singular \n")
        }
    }
    
    ### ** Compute derivative
    if(is.null(object$dVcov)){
        dVcov.dtheta  <- dVcov2(object, ...)
    }else{
        dVcov.dtheta <- object$dVcov
    }
    vcov.param <- attr(dVcov.dtheta, "vcov.param")
    attr(dVcov.dtheta, "vcov.param") <- NULL
    keep.param <- dimnames(dVcov.dtheta)[[3]]
    
    ### ** Compute degrees of freedom
    calcDF <- function(M.C){ # M.C <- C
        C.vcov.C <- rowSums(M.C %*% vcov.param * M.C)
    
        C.dVcov.C <- sapply(keep.param, function(x){
            rowSums(M.C %*% dVcov.dtheta[,,x] * M.C)
        })
        
        numerator <- 2 *(C.vcov.C)^2
        denom <- rowSums(C.dVcov.C %*% vcov.param[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
        df <- numerator/denom
        return(df)
    }

    df.Wald  <- calcDF(C)
    
    ### ** Anova
    ## F test
    ## Fstat <- t(C %*% p) %*% i.C.vcov.C %*% (C %*% p) / q

    ## df
    svd.tempo <- eigen(solve(C %*% vcov.param %*% t(C)))
    D.svd <- diag(svd.tempo$values)
    P.svd <- svd.tempo$vectors

    C.anova <- sqrt(D.svd) %*% t(P.svd) %*% C
    ## Fstat - crossprod(C.anova %*% p)/q
    df.F <- calcDF(C.anova)

    ## ** export
    return(list(Wald = df.Wald,
                F = df.F))
    
}



##----------------------------------------------------------------------
### dfVariance.R ends here
