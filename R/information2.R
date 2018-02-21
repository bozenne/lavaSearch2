### information2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (14:17) 
## Version: 
## Last-Updated: feb 21 2018 (18:05) 
##           By: Brice Ozenne
##     Update #: 50
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .information2
#' @title Compute the Expected Information Matrix From the Conditional Moments
#' @description Compute the expected information matrix from the conditional moments
#' when there is no missing value.
#'
#' @param n.corrected [numeric >0] the number of independent epsilon.
#' May not match the number of observations when a small sample correction is used.
#'
#' @details \code{.information2} will perform the computation individually when the
#' argument \code{index.Omega} is null.
#' 
#' @keywords internal
.information2 <- function(dmu, dOmega,
                          Omega, OmegaM1, n.corrected, leverage, index.Omega, n.cluster,
                          grid.meanparam, n.grid.meanparam,
                          grid.varparam, n.grid.varparam,
                          name.param, name.meanparam, name.varparam,
                          param2index, n.param){

    ### ** Prepare
    test.global <- is.null(index.Omega)
    
    Info <- matrix(0, nrow = n.param, ncol = n.param,
                   dimnames = list(name.param,name.param))

    ### ** Global
    if(test.global){
        ### *** Information relative to the mean parameters
        for(iG in 1:n.grid.meanparam){ # iG <- 1
            iName1 <- grid.meanparam[iG,1]
            iName2 <- grid.meanparam[iG,2]
            iP1 <- param2index[iName1]
            iP2 <- param2index[iName2]

            Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iName1]] %*% OmegaM1 * dmu[[iName2]])
        }

        ### *** Information realtive to the variance parameters
        for(iG in 1:n.grid.varparam){ # iG <- 1
            iName1 <- grid.varparam[iG,1]
            iName2 <- grid.varparam[iG,2]
            iP1 <- param2index[iName1]
            iP2 <- param2index[iName2]

            iDiag <- diag(OmegaM1 %*% dOmega[[iName1]] %*% OmegaM1 %*% dOmega[[iName2]])
            Info[iP1,iP2] <- Info[iP1,iP2] + 1/2*sum(iDiag*n.corrected)
        }

    }
    
    ### ** Individual specific
    if(!test.global){
        for(iC in 1:n.cluster){ # iC <- 1
            iIndex <- index.Omega[[iC]]

            ### *** Information relative to the mean parameters
            for(iG in 1:n.grid.meanparam){ # iG <- 1
                iName1 <- grid.meanparam[iG,1]
                iName2 <- grid.meanparam[iG,2]
                iP1 <- param2index[iName1]
                iP2 <- param2index[iName2]

                Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iName1]][iC,iIndex] %*% OmegaM1[iIndex,iIndex,drop=FALSE] * dmu[[iName2]][iC,iIndex])            
            }

            ### *** Information realtive to the variance parameters
            for(iG in 1:n.grid.varparam){ # iG <- 1
                iName1 <- grid.varparam[iG,1]
                iName2 <- grid.varparam[iG,2]
                iP1 <- param2index[iName1]
                iP2 <- param2index[iName2]

                iDiag <- diag(OmegaM1[iIndex,iIndex,drop=FALSE] %*% dOmega[[iName1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[iIndex,iIndex,drop=FALSE] %*% dOmega[[iName2]][iIndex,iIndex,drop=FALSE])
                Info[iP1,iP2] <- Info[iP1,iP2] + 1/2 * sum(iDiag * (1 - leverage[iC,iIndex]))            
            }
        }        
    }


    ### ** Make Info a symmetric matrix
    Info <- symmetrize(Info, update.upper = NULL)
    
    ### ** export
    return(Info)
}

##----------------------------------------------------------------------
### information2.R ends here
