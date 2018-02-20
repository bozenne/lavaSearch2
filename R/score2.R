### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: feb 19 2018 (18:10) 
##           By: Brice Ozenne
##     Update #: 2219
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
                    name.param, name.meanparam, name.varparam, n.cluster,
                    indiv){

    
    out.score <- matrix(0, nrow = n.cluster, ncol = length(name.param),
                        dimnames = list(NULL,name.param))
    epsilon.OmegaM1 <- epsilon %*% OmegaM1

    ## ** Compute score relative to the mean coefficients
    for(iP in name.meanparam){ # iP <- 1
        out.score[,iP] <- out.score[,iP] + rowSums(dmu[[iP]] * epsilon.OmegaM1)            
    }
    ## ** Compute score relative to the variance-covariance coefficients
    for(iP in name.varparam){ # iP <- 1
        term2 <- - 1/2 * tr(OmegaM1 %*% dOmega[[iP]])            
        term3 <- 1/2 * rowSums(epsilon.OmegaM1 %*% dOmega[[iP]] * epsilon.OmegaM1)
        out.score[,iP] <- out.score[,iP] + as.double(term2) + term3 
    }        

    ### ** export
    if(indiv==FALSE){
        out.score <- colSums(out.score)
    }
    return(out.score)
}

## * .score2Indiv
#' @title Compute the Corrected Score.
#' @description Compute the corrected score in presence of missing values.
#'
#' @param n.cluster [interger >0] the number of observations.
#' 
#' @keywords internal
.score2Indiv <- function(dmu.dtheta, dOmega.dtheta, epsilon,
                    Omega, ls.indexOmega,
                    indiv, 
                    name.param, n.param, n.cluster){

    clusterSpecific <- !is.null(ls.indexOmega)
    name.meanparam <- names(dmu.dtheta)
    name.vcovparam <- names(dOmega.dtheta)
    out.score <- matrix(0, nrow = n.cluster, ncol = n.param,
                        dimnames = list(NULL,name.param))

### ** Individual specific Omega (e.g. presence of missing values)
    if(clusterSpecific){
        for(iC in 1:n.cluster){
            iOmega.tempo <- chol2inv(chol(Omega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]))
            epsilon.iOmega.tempo <- iOmega.tempo %*% cbind(epsilon[iC,ls.indexOmega[[iC]]])

            ## *** Compute score relative to the mean coefficients
            for(iP in name.meanparam){ # iP <- name.meanparam[1]
                out.score[iC,iP] <- out.score[iC,iP] + dmu.dtheta[[iP]][iC,ls.indexOmega[[iC]]] %*% epsilon.iOmega.tempo
            }

            ## *** Compute score relative to the variance-covariance coefficients
            for(iP in name.vcovparam){ # iP <- name.vcovparam[1]
                dOmega.dtheta.tempo <-  dOmega.dtheta[[iP]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                                                            
                term2 <- - 1/2 * tr(iOmega.tempo %*% dOmega.dtheta.tempo)
                term3 <- 1/2 * sum(epsilon.iOmega.tempo * dOmega.dtheta.tempo %*% epsilon.iOmega.tempo)
                out.score[iC,iP] <- out.score[iC,iP] + as.double(term2) + term3 
            }
        }
    }
            
    ### ** Same for all individuals
    if(clusterSpecific == FALSE){
        iOmega <- chol2inv(chol(Omega))
        epsilon.iOmega <- epsilon %*% iOmega

        ## *** Compute score relative to the mean coefficients
        for(iP in name.meanparam){ # iP <- 1
            out.score[,iP] <- out.score[,iP] + rowSums(dmu.dtheta[[iP]] * epsilon.iOmega)            
        }
        ## *** Compute score relative to the variance-covariance coefficients
        for(iP in name.vcovparam){ # iP <- 1
            term2 <- - 1/2 * tr(iOmega %*% dOmega.dtheta[[iP]])            
            term3 <- 1/2 * rowSums(epsilon.iOmega %*% dOmega.dtheta[[iP]] * epsilon.iOmega)
            out.score[,iP] <- out.score[,iP] + as.double(term2) + term3 
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
