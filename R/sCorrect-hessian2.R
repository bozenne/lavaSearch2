### sCorrect-hessian2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (14:09) 
## Version: 
## Last-Updated: jan 10 2020 (13:33) 
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

## * .hessian2
#' @title Compute the Hessian Matrix From the Conditional Moments
#' @description Compute the Hessian matrix from the conditional moments.
#' @name hessian2-internal
#' 
#' @details \code{calc_hessian} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.hessian2 <- function(dmu, dOmega, d2mu, d2Omega, epsilon, OmegaM1,
                      missing.pattern, unique.pattern, name.pattern,
                      grid.mean, grid.var, grid.hybrid, name.param,
                      leverage, n.cluster){
    if(lava.options()$debug){cat(".hessian2\n")}

    ## ** Prepare
    n.grid.mean <- NROW(grid.mean)
    n.grid.var <- NROW(grid.var)
    n.grid.hybrid <- NROW(grid.hybrid)
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)

    hessian <- array(0, dim = c(n.param, n.param, n.cluster),
                     dimnames = list(name.param,name.param,NULL))
    if(length(dmu)>0){
        index.mean <- 1:n.grid.mean
    }else{
        index.mean <- NULL
    }
    if(length(dOmega)>0){
        index.var <- 1:n.grid.var
    }else{
        index.var <- NULL
    }
    if(length(dmu)>0 && length(dOmega)>0){
        index.hybrid <- 1:n.grid.hybrid
    }else{
        index.hybrid <- NULL
    } 

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        iEpsilon <- epsilon[iIndex,iY,drop=FALSE]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)
        iLeverage <- leverage[iIndex,iY,drop=FALSE]

        ## *** second derivative relative to the mean parameters
        for(iG in index.mean){ # iG <- 2
            iP1 <- grid.mean[iG,1]
            iP2 <- grid.mean[iG,2]
            if(grid.mean[iG,"d2.12"]){
                term1 <- rowSums((id2mu[[iP1]][[iP2]] %*% iOmegaM1) * iEpsilon)
            }else if(grid.mean[iG,"d2.21"]){
                term1 <- rowSums((id2mu[[iP2]][[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1 <- 0
            }
            term2 <- -rowSums((idmu[[iP1]] %*% iOmegaM1) * idmu[[iP2]])
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] + term1 + term2
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
        ## *** second derivative relative to the variance parameters
        for(iG in index.var){ # iG <- 2
            iP1 <- grid.var[iG,1]
            iP2 <- grid.var[iG,2]

            term1a <- - diag(iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1 %*% idOmega[[iP2]])
            term2 <- - rowSums((iEpsilon %*% iOmegaM1 %*% idOmega[[iP2]] %*% iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1) * iEpsilon)
            if(grid.var[iG,"d2.12"]){
                term1b <- diag(iOmegaM1 %*% id2Omega[[iP1]][[iP2]])
                term3 <- 1/2 * rowSums((iEpsilon %*% iOmegaM1 %*% id2Omega[[iP1]][[iP2]] %*% iOmegaM1) * iEpsilon)
            }else if(grid.var[iG,"d2.21"]){
                term1b <- diag(iOmegaM1 %*% id2Omega[[iP2]][[iP1]])
                term3 <- 1/2 * rowSums((iEpsilon %*% iOmegaM1 %*% id2Omega[[iP2]][[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1b <- 0
                term3 <- 0
            }
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] - 1/2 * rowSums( sweep(1-iLeverage, FUN = "*", STATS = term1a + term1b, MARGIN = 2) ) + term2 + term3
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
        
        ## *** second derivative relative to the mean and variance parameters
        for(iG in index.hybrid){ # iG <- 1
            iP1 <- grid.hybrid[iG,1]
            iP2 <- grid.hybrid[iG,2]

            if(!is.null(idmu[[iP1]]) && !is.null(idOmega[[iP2]])){
                term1 <- - rowSums((idmu[[iP1]] %*% iOmegaM1 %*% idOmega[[iP2]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1 <- 0
            }
            if(!is.null(idmu[[iP2]]) && !is.null(idOmega[[iP1]])){
                term2 <- - rowSums((idmu[[iP2]] %*% iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term2 <- 0
            }
            
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] + term1 + term2
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
    }

    ## ** export
    return(hessian)
}

## * .subsetList
.subsetList <- function(object, indexRow, indexCol){
    if(length(object)==0){
        return(object)
    }else{    
        return(lapply(object, FUN = function(iL){iL[indexRow,indexCol,drop=FALSE]}))
    }
}
## * .subsetList2
.subsetList2 <- function(object, indexRow, indexCol){
    if(length(object)==0){
        return(object)
    }else{    
        return(lapply(object, FUN = .subsetList, indexRow = indexRow, indexCol = indexCol))
    }
}


######################################################################
### sCorrect-hessian2.R ends here
