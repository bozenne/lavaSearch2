### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: feb 21 2018 (18:10) 
##           By: Brice Ozenne
##     Update #: 2233
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * .score2
#' @title Compute the Corrected Score.
#' @description Compute the corrected score when there is no missing value.
#'
#' @param n.cluster [interger >0] the number of observations.
#' 
#' @keywords internal
.score2 <- function(epsilon, Omega, OmegaM1, dmu, dOmega,                    
                    name.param, name.meanparam, name.varparam,
                    index.Omega, n.cluster, indiv){

    ## ** Prepare
    test.global <- is.null(index.Omega)
    out.score <- matrix(0, nrow = n.cluster, ncol = length(name.param),
                        dimnames = list(NULL,name.param))
    
    
            
    ### ** global
    if(test.global){
        epsilon.OmegaM1 <- epsilon %*% OmegaM1

        ## *** Compute score relative to the mean coefficients
        for(iP in name.meanparam){ # iP <- 1
            out.score[,iP] <- out.score[,iP] + rowSums(dmu[[iP]] * epsilon.OmegaM1)
        }
        
        ## *** Compute score relative to the variance-covariance coefficients
        for(iP in name.varparam){ # iP <- 1
            term2 <- - 1/2 * tr(OmegaM1 %*% dOmega[[iP]])            
            term3 <- 1/2 * rowSums(epsilon.OmegaM1 %*% dOmega[[iP]] * epsilon.OmegaM1)
            out.score[,iP] <- out.score[,iP] + as.double(term2) + term3
        }        
    }

    ### ** individual specific
    if(!test.global){
        for(iC in 1:n.cluster){
            iIndex <- index.Omega[[iC]]
            
            iOmegaM1 <- chol2inv(chol(Omega[iIndex,iIndex,drop=FALSE]))
            iEpsilon.OmegaM1 <- iOmegaM1 %*% cbind(epsilon[iC,iIndex])


            ## *** Compute score relative to the mean coefficients
            for(iP in name.meanparam){ # iP <- name.meanparam[1]
                out.score[iC,iP] <- out.score[iC,iP] + dmu[[iP]][iC,iIndex] %*% iEpsilon.OmegaM1
            }

            ## *** Compute score relative to the variance-covariance coefficients
            for(iP in name.varparam){ # iP <- name.varparam[1]
                term2 <- - 1/2 * tr(iOmegaM1 %*% dOmega[[iP]][iIndex,iIndex,drop=FALSE])
                term3 <- 1/2 * sum(iEpsilon.OmegaM1 * dOmega[[iP]][iIndex,iIndex,drop=FALSE] %*% iEpsilon.OmegaM1)
                out.score[iC,iP] <- out.score[iC,iP] + as.double(term2) + term3 
            }
        }
    }

    ### ** export
    if(indiv==FALSE){
        out.score <- colSums(out.score)
    }
    return(out.score)
}


#----------------------------------------------------------------------
### score2.R ends her
