### biasCoxSnell.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (10:20) 
## Version: 
## Last-Updated: jan  9 2020 (18:01) 
##           By: Brice Ozenne
##     Update #: 49
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .init_sscCoxSnell
.init_sscCoxSnell <- function(object,...){
    return(object,...)
}

## * .sscCoxSnell
.sscCoxSnell <- function(object){

    ## ** compute JJK
    JJK <- .calcJJK(object)


    
    names(object$sCorrect$dmoment)
    names(object$sCorrect$skeleton)
    object$sCorrect$skeleton$grid.dmoment
    
    param <- coef2(e)
    name.param <- names(param)
    n.param <- length(name.param)
    data <- model.frame(e)
    n <- nobs(e)

    ## *** identify missing values
    names(object)
    object$cluster
    
    ## *** pre-compute skeleton for conditional moments
    if(is.null(object$sCorrect$conditionalMoment)){
        object$sCorrect$conditionalMoment <- conditionalMoment(e, data = data,
                                                               first.order = TRUE, second.order = TRUE,
                                                               name.endogenous = endogenous(e),
                                                               name.latent = latent(e), usefit = FALSE)
    }

    
    ## ** update conditional moments
    moment.object <- conditionalMoment(e, data = data, param = param,
                                       first.order = TRUE, second.order = TRUE,
                                       name.endogenous = endogenous(e),
                                       name.latent = latent(e), usefit = TRUE)

    Omega <- moment.object$Omega
    dmu <- moment.object$dmu
    dOmega <- moment.object$dOmega
    d2mu <- moment.object$d2mu
    d2Omega <- moment.object$d2Omega

    ## ** update variance-covariance matrix of the parameters
    if(is.null(vcov)){
        grid.information <- .gridInformation(param = name.param,
                                             param.mean = attr(param, "mean.coef"),
                                             param.var = attr(param, "var.coef"))
        
        iVcov <- .vcov2(dmu = dmu,
                        dOmega = dOmega,
                        Omega = Omega,
                        n.corrected = n,
                        leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                        grid.meanparam = grid.information$mean,
                        n.grid.meanparam = grid.information$n.mean,
                        grid.varparam = grid.information$var,
                        n.grid.varparam = grid.information$n.var,
                        name.param = name.param,
                        n.param = n.param,
                        attr.info = TRUE)
        iInfo <- .information2(dmu = dmu,
                               dOmega = dOmega,
                               Omega = Omega,
                               n.corrected = n.corrected,
                               leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                               grid.meanparam = grid.meanparam,
                               n.grid.meanparam = n.grid.meanparam,
                               grid.varparam = grid.varparam,
                               n.grid.varparam = n.grid.varparam,
                               name.param = name.param,
                               n.param = n.param)
        iVcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
        if(inherits(iVcov.param, "try-error")){
            iVcov.param <- solve(iInfo)
        }

    }
    

    ## 2*J(\theta3;\theta2,\theta1) + K(\theta3,\theta2,\theta1)    

    moment.object$value
}


## * .calcJJK
.calcJJK <- function(object){

    ## ** extract information
    Omega <- object$sCorrect$moment$Omega
    dmu <- object$sCorrect$dmoment$dmu
    d2mu <- object$sCorrect$d2moment$d2mu
    dOmega <- object$sCorrect$dmoment$dOmega
    d2Omega <- object$sCorrect$d2moment$d2Omega

    missing.pattern <- object$sCorrect$missing$pattern
    name.pattern <- object$sCorrect$missing$name.pattern
    unique.pattern <- object$sCorrect$missing$unique.pattern
    n.pattern <- length(name.pattern)
    OmegaM1 <- object$sCorrect$moment$Omega.missing.pattern
    
    name.param <- names(object$sCorrect$param)
    n.param <- length(name.param)

    grid.2meanD1.1varD1 <- object$sCorrect$skeleton$grid.2meanD1.1varD1
    grid.2meanD2.1meanD1 <- object$sCorrect$skeleton$grid.2meanD2.1meanD1
    grid.2varD2.1varD1 <- object$sCorrect$skeleton$grid.2varD2.1varD1
    n.grid.2meanD1.1varD1 <- NROW(grid.2meanD1.1varD1)
    n.grid.2meanD2.1meanD1 <- NROW(grid.2meanD2.1meanD1)
    n.grid.2varD2.1varD1 <- NROW(grid.2varD2.1varD1)

    ## ** prepare output    
    JJK <-  array(0, dim = c(n.param,n.param,n.param),
                  dimnames = list(name.param,name.param,name.param))

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)
        
        ## *** 1 second derivative and 1 first derivative regarding the variance
        if(n.grid.2varD2.1varD1>0){
            for(iGrid in 1:n.grid.2varD2.1varD1){ # iGrid <- 1
                iName1 <- grid.2varD2.1varD1[iGrid,"X"]
                iName2 <- grid.2varD2.1varD1[iGrid,"Y"]
                iName3 <- grid.2varD2.1varD1[iGrid,"Z"]

                ## term 1
                iDiag <- diag(iOmegaM1 %*% id2Omega[[iName2]][[iName3]] %*% iOmegaM1 %*% idOmega[[iName1]])
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + sum(iDiag * iN.corrected)

                ## term 2
                iDiag <- diag(iOmegaM1 %*% id2Omega[[iName2]][[iName3]] %*% iOmegaM1 %*% idOmega[[iName1]])
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + sum(iDiag * iN.corrected)

                ## term 3
                iDiag <- - diag(iOmegaM1 %*% id2Omega[[iName2]][[iName3]] %*% iOmegaM1 %*% idOmega[[iName1]])
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + sum(iDiag * iN.corrected)

                ## symmetrize
                dInfo[iName2,iName1,iName3] <- dInfo[iName1,iName2,iName3]
            }
        }
        
        ## *** 2 first derivative regarding the mean and one regarding the variance
        if(n.grid.2meanD1.1varD1>0){
            for(iGrid in 1:n.grid.2meanD1.1varD1){ # iGrid <- 1
                iName1 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"Y"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Z"]
            
                ## term 4
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] - sum(idmu[[iName1]] %*% iOmegaM1 %*% idOmega[[iName3]] %*% iOmegaM1 * idmu[[iName2]])
                dInfo[iName2,iName1,iName3] <- dInfo[iName1,iName2,iName3] ## symmetrize
            }
        }
        
        ## *** 1 second derivative and 1 first derivative regarding the mean
        if(n.grid.2meanD2.1meanD1>0){
            for(iGrid in 1:n.grid.2meanD2.1meanD1){ # iGrid <- 1
                iName1 <- grid.2meanD2.1meanD1[iGrid,"X"]
                iName2 <- grid.2meanD2.1meanD1[iGrid,"Y"]
                iName3 <- grid.2meanD2.1meanD1[iGrid,"Z"]

                ## term 5
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + sum(id2mu[[iName2]][[iName3]] %*% iOmegaM1 * idmu[[iName1]])
                dInfo[iName2,iName1,iName3] <- dInfo[iName2,iName1,iName3] ## symmetrize

                ## tern 6
                dInfo[iName1,iName2,iName3] <- dInfo[iName1,iName2,iName3] + sum(idmu[[iName1]] %*% iOmegaM1 * id2mu[[iName2]][[iName3]])
                dInfo[iName2,iName1,iName3] <- dInfo[iName1,iName2,iName3] ## symmetrize
            }
        }
    }

    return(JJK)
}

######################################################################
### biasCoxSnell.R ends here
