### sCorrect-dInformation2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (14:09) 
## Version: 
## Last-Updated: dec 12 2019 (09:07) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .dInformation2
#' @title Compute the First Derivative of the Expected Information Matrix
#' @description Compute the first derivative of the expected information matrix.
#' @name .dinformation2-internal
#' 
#' @details \code{calc_dinformation} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.dInformation2 <- function(dmu, dOmega, d2mu, d2Omega, OmegaM1,
                           missing.pattern, unique.pattern, name.pattern,
                           grid.meanparam, grid.varparam, grid.hybrid, name.param, name.3deriv,
                           leverage, n.cluster){

    if(TRUE){cat(".dInformation2\n")}

    ## ** Prepare
    n.grid.meanparam <- NROW(grid.meanparam)
    n.grid.varparam <- NROW(grid.varparam)
    n.grid.hybridparam <- NROW(grid.hybridparam)
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)
    n.3deriv <- length(name.3deriv)
        
    n.param <- length(name.param)
    index.deriv <- match(name.3deriv, name.param)

    dInfo <-  array(0,
                    dim = c(n.param, n.param, n.3deriv),
                    dimnames = list(name.param, name.param, name.3deriv))

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iN.corrected <- n.cluster - colSums(leverage[iIndex,iY,drop=FALSE])

        for(iDeriv in index.deriv){ # iDeriv <- 4
            for(iP1 in 1:n.param){ # iP1 <- 1
                for(iP2 in iP1:n.param){ # iP2 <- 1
                
                    iNameD <- name.param[iDeriv]
                    iName1 <- name.param[iP1]
                    iName2 <- name.param[iP2]
                    ## cat(iNameD," ",iName1,"",iName2,"\n")
                
                    ## *** identify relevant terms
                    test.Omega1 <- !is.null(dOmega[[iNameD]]) && !is.null(dOmega[[iName1]]) && !is.null(dOmega[[iName2]])
                    test.Omega2a <- !is.null(d2Omega[[iNameD]][[iName1]]) && !is.null(dOmega[[iName2]])
                    test.Omega2b <- !is.null(d2Omega[[iName1]][[iNameD]]) && !is.null(dOmega[[iName2]])
                    test.Omega3a <- !is.null(d2Omega[[iNameD]][[iName2]]) && !is.null(dOmega[[iName1]])
                    test.Omega3b <- !is.null(d2Omega[[iName2]][[iNameD]]) && !is.null(dOmega[[iName1]])
                
                    test.mu1a <- !is.null(d2mu[[iNameD]][[iName1]]) && !is.null(dmu[[iName2]])
                    test.mu1b <- !is.null(d2mu[[iName1]][[iNameD]]) && !is.null(dmu[[iName2]])
                    test.mu2a <- !is.null(d2mu[[iNameD]][[iName2]]) && !is.null(dmu[[iName1]])
                    test.mu2b <- !is.null(d2mu[[iName2]][[iNameD]]) && !is.null(dmu[[iName1]])
                    test.mu3 <- !is.null(dOmega[[iNameD]]) && !is.null(dmu[[iName1]]) && !is.null(dmu[[iName2]])

                    if((test.Omega1 + test.Omega2a + test.Omega2b + test.Omega3a + test.Omega3b + test.mu1a + test.mu1b + test.mu2a + test.mu2b + test.mu3) == 0){
                        next
                    }

                    ## *** extract quantities for computations 
                    if(test.mu1a){
                        d2mu.D1 <- d2mu[[iNameD]][[iName1]][iIndex,iY,drop=FALSE]
                    }else if(test.mu1b){
                        d2mu.D1 <- d2mu[[iName1]][[iNameD]][iIndex,iY,drop=FALSE]
                    }
                    if(test.mu2a){
                        d2mu.D2 <- d2mu[[iNameD]][[iName2]][iIndex,iY,drop=FALSE]
                    }else if(test.mu2b){
                        d2mu.D2 <- d2mu[[iName2]][[iNameD]][iIndex,iY,drop=FALSE]
                    }
                    if(test.Omega2a){
                        d2Omega.D1 <- d2Omega[[iNameD]][[iName1]][iY,iY,drop=FALSE]
                    }else if(test.Omega2b){
                        d2Omega.D1 <- d2Omega[[iName1]][[iNameD]][iY,iY,drop=FALSE]
                    }
                    if(test.Omega3a){
                        d2Omega.D2 <- d2Omega[[iNameD]][[iName2]][iY,iY,drop=FALSE]
                    }else{
                        d2Omega.D2 <- d2Omega[[iName2]][[iNameD]][iY,iY,drop=FALSE]
                    }
                
                    ## *** pre-compute 
                    if(!is.null(dOmega[[iNameD]])){
                        OmegaM1.dOmega.D <- iOmegaM1 %*% dOmega[[iNameD]][iY,iY,drop=FALSE]
                    }
                    if(!is.null(dOmega[[iName1]])){
                        OmegaM1.dOmega.1 <- iOmegaM1 %*% dOmega[[iName1]][iY,iY,drop=FALSE]
                    }
                    if(!is.null(dOmega[[iName2]])){
                        OmegaM1.dOmega.2 <- iOmegaM1 %*% dOmega[[iName2]][iY,iY,drop=FALSE]
                    }                    

                    ## *** evaluate contributions to dInformation
                    if(test.Omega1){
                        iDiag1 <- diag(OmegaM1.dOmega.D %*% OmegaM1.dOmega.1 %*% OmegaM1.dOmega.2)
                        iDiag2 <- diag(OmegaM1.dOmega.1 %*% OmegaM1.dOmega.D %*% OmegaM1.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iN.corrected + iDiag2 * iN.corrected)
                    }

                    if(test.Omega2a || test.Omega2b){
                        iDiag <- diag(OmegaM1 %*% d2Omega.D1 %*% OmegaM1.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.corrected)
                    }

                    if(test.Omega3a || test.Omega3b){
                        iDiag <- diag(OmegaM1.dOmega.1 %*% OmegaM1 %*% d2Omega.D2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.corrected)
                    }

                    if(test.mu1a || test.mu1b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% OmegaM1 * dmu[[iName2]][iIndex,iY,drop=FALSE])
                    }

                    if(test.mu2a || test.mu2b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu[[iName1]][iIndex,iY,drop=FALSE] %*% OmegaM1 * d2mu.D2)
                    }
                    
                    if(test.mu3){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu[[iName1]][iIndex,iY,drop=FALSE] %*% OmegaM1.dOmega.D %*% OmegaM1 * dmu[[iName2]][iIndex,iY,drop=FALSE])
                    }
                }
            }
            
        }
    }

    ## ** Symmetrize
    dInfo[,,iNameD] <- symmetrize(dInfo[,,iNameD], update.upper = NULL)

    ### ** export
    return(dInfo)
}


######################################################################
### sCorrect-dInformation2.R ends here
