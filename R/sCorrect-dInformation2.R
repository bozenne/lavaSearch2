### sCorrect-dInformation2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (14:09) 
## Version: 
## Last-Updated: jan  7 2020 (13:46) 
##           By: Brice Ozenne
##     Update #: 62
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
                           grid.dInformation, name.param, name.param.dInformation,
                           leverage, n.cluster){

    if(lava.options()$debug){cat(".dInformation2\n")}

    ## ** Prepare
    n.grid <- NROW(grid.dInformation)
    n.param <- length(name.param)
    n.param.dInformation <- length(name.param.dInformation)
    n.pattern <- length(name.pattern)

    dInfo <-  array(0,
                    dim = c(n.param, n.param, n.param.dInformation),
                    dimnames = list(name.param, name.param, name.param.dInformation))

    index.duplicated <- which(grid.dInformation$duplicated)
    index.Nduplicated <- setdiff(1:n.grid, index.duplicated)
    ## grid.dInformation[index.Nduplicated,,drop=FALSE]
    
    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iN.corrected <- n.cluster - colSums(leverage[iIndex,iY,drop=FALSE])
        for(iGrid in index.Nduplicated){ # iGrid <- 1
            iName1 <- grid.dInformation[iGrid,"X"]
            iName2 <- grid.dInformation[iGrid,"Y"]
            iNameD <- grid.dInformation[iGrid,"Z"]
            ## cat("* ", iNameD," ",iName1,"",iName2,"\n")

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
                iOmegaM1.dOmega.D <- iOmegaM1 %*% dOmega[[iNameD]][iY,iY,drop=FALSE]
            }
            if(!is.null(dOmega[[iName1]])){
                iOmegaM1.dOmega.1 <- iOmegaM1 %*% dOmega[[iName1]][iY,iY,drop=FALSE]
            }
            if(!is.null(dOmega[[iName2]])){
                iOmegaM1.dOmega.2 <- iOmegaM1 %*% dOmega[[iName2]][iY,iY,drop=FALSE]
            }                    

            ## *** evaluate contributions to dInformation
            ## if(iP==2 && (iName1==iName2)&& (iName2==iNameD) && (iName1=="2")){browser()}
            if(test.Omega1){                
                iDiag1 <- diag(iOmegaM1.dOmega.D %*% iOmegaM1.dOmega.1 %*% iOmegaM1.dOmega.2)
                iDiag2 <- diag(iOmegaM1.dOmega.1 %*% iOmegaM1.dOmega.D %*% iOmegaM1.dOmega.2)
                dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iN.corrected + iDiag2 * iN.corrected)
            }

            if(test.Omega2a || test.Omega2b){
                iDiag <- diag(iOmegaM1 %*% d2Omega.D1 %*% iOmegaM1.dOmega.2)
                dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.corrected)
            }

            if(test.Omega3a || test.Omega3b){
                iDiag <- diag(iOmegaM1.dOmega.1 %*% iOmegaM1 %*% d2Omega.D2)
                dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.corrected)
                
            }

            if(test.mu1a || test.mu1b){
                dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% iOmegaM1 * dmu[[iName2]][iIndex,iY,drop=FALSE])
            }

            if(test.mu2a || test.mu2b){
                dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu[[iName1]][iIndex,iY,drop=FALSE] %*% iOmegaM1 * d2mu.D2)
            }
                    
            if(test.mu3){
                dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu[[iName1]][iIndex,iY,drop=FALSE] %*% iOmegaM1.dOmega.D %*% iOmegaM1 * dmu[[iName2]][iIndex,iY,drop=FALSE])
            }

        }
    }

    ## ** symmetrize
    if(length(index.duplicated)>0){
        for(iGrid in index.duplicated){ ## iGrid <- index.duplicated[1]


            iName1 <- grid.dInformation[iGrid,"X"]
            iName2 <- grid.dInformation[iGrid,"Y"]
            iNameD <- grid.dInformation[iGrid,"Z"]

            iName1.ref <- grid.dInformation[grid.dInformation[iGrid,"reference"],"X"]
            iName2.ref <- grid.dInformation[grid.dInformation[iGrid,"reference"],"Y"]
            iNameD.ref <- grid.dInformation[grid.dInformation[iGrid,"reference"],"Z"]
            
            dInfo[iName1,iName2,iNameD] <- dInfo[iName1.ref,iName2.ref,iNameD.ref]
        }
    }

    ### ** export
    return(dInfo)
}


######################################################################
### sCorrect-dInformation2.R ends here
