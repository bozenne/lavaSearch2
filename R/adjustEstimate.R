### adjustEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: feb 20 2018 (17:43) 
##           By: Brice Ozenne
##     Update #: 305
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * adjustEstimate
#' @title Compute bias corrected quantities.
#' @description Compute bias corrected quantities when there is no missing value.
#' 
#' @keywords internal
adjustEstimate <- function(epsilon, Omega, dmu, dOmega, n.cluster,
                           name.param, name.endogenous, name.meanparam, name.varparam,
                           n.endogenous.cluster, index.Omega,
                           adjust.Omega, adjust.n, tol, n.iter){

    ## ** Prepare
    name.hybridparam <- intersect(name.meanparam, name.varparam)

    n.param <- length(name.param)
    n.meanparam <- length(name.meanparam)
    n.varparam <- length(name.varparam)
    n.hybridparam <- length(name.hybridparam)
    n.endogenous <- length(name.endogenous)

    grid.meanparam <- .combination(name.meanparam, name.meanparam)
    n.grid.meanparam <- NROW(grid.meanparam)
    grid.varparam <- .combination(name.varparam, name.varparam)
    n.grid.varparam <- NROW(grid.varparam)

    param2index <- setNames(1:n.param, name.param)


    if(is.null(index.Omega)){
        n.endogenous.cluster <- rep(n.endogenous, n.cluster)
        n.corrected <- rep(n.cluster, n.endogenous)        
        epsilon <- sapply(1:n.cluster, function(iC){list(epsilon[iC,])})
    }else{
        n.corrected <- NULL        
    }

    ls.Psi <- vector(mode = "list", length = n.cluster)
    leverage <- vector(mode = "list", length = n.cluster)
    ls.dmu <- vector(mode = "list", length = n.cluster)
    for(iC in 1:n.cluster){ # iC <- 1
        ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = n.endogenous.cluster[[iC]],
                               dimnames = list(name.param, name.endogenous[index.Omega[[iC]]]))
        if(is.null(index.Omega)){
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu,function(x){x[iC,]}))
        }else{
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu,function(x){x[[iC]]}))
        }
        leverage[[iC]] <- rep(0,n.endogenous.cluster[[iC]])
    }
    
    ## ** Initialisation (i.e. first iteration without correction)
    Omega.adj <- Omega
    OmegaM1.adj <- chol2inv(chol(Omega.adj))

    iInfo <- .information2(dmu = dmu,
                           dOmega = dOmega,
                           Omega = Omega.adj,
                           OmegaM1 = OmegaM1.adj,
                           n.corrected = n.corrected,
                           leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                           grid.meanparam = grid.meanparam,
                           n.grid.meanparam = n.grid.meanparam,
                           grid.varparam = grid.varparam,
                           n.grid.varparam = n.grid.varparam,
                           name.param = name.param,
                           name.meanparam = name.meanparam,
                           name.varparam = name.varparam,
                           param2index = param2index, n.param = n.param)
    iVcov.param <- solve(iInfo)
    epsilon.adj <- epsilon
    
    ## ** Loop
    if(adjust.Omega || adjust.n){
        iIter <- 0
        iTol <- Inf
    }else{
        iIter <- Inf
        iTol <- -Inf
    }
    
    while(iIter < n.iter & iTol > tol){
        Omega_save <- Omega.adj

        ## *** Step (ii-iii): compute individual bias and expected bias
        Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                      dimnames = list(name.endogenous, name.endogenous))
        for(iC in 1:n.cluster){
            ## individual bias        
            ls.Psi[[iC]] <- t(ls.dmu[[iC]])  %*% iVcov.param %*% ls.dmu[[iC]]

            ## average bias            
            if(is.null(index.Omega)){
                Psi <- Psi + ls.Psi[[iC]]/n.cluster                
            }else{
                Psi[index.Omega[[iC]],index.Omega[[iC]]] <- Psi[index.Omega[[iC]],index.Omega[[iC]]] + ls.Psi[[iC]]/n.cluster
            }
        }

        ## *** Step (iv): correct residual covariance matrix
        if(adjust.Omega){
            ## corrected residual covariance variance
            Omega.adj <- Omega + Psi
            OmegaM1.adj <- chol2inv(chol(Omega.adj))
        }

        ## *** Step (v): corrected sample size
        if(adjust.n){
            Omega.adj.chol <- chol(Omega.adj)
            if(is.null(index.Omega)){
                iOmega.adj <- Omega.adj
                iOmegaM1.adj <- OmegaM1.adj
                iOmega.adj.chol <- Omega.adj.chol
                iIndex.Omega <- 1:n.endogenous
            }
            
            for(iC in 1:n.cluster){                 # iC <- 1
                if(!is.null(index.Omega)){
                    iOmega.adj <- Omega.adj[index.Omega[[iC]],index.Omega[[iC]]]
                    iOmegaM1.adj <- OmegaM1.adj[index.Omega[[iC]],index.Omega[[iC]]]
                    iOmega.adj.chol <- Omega.adj.chol[index.Omega[[iC]],index.Omega[[iC]]]
                    iIndex.Omega <- index.Omega[[iC]]
                }
                browser()
                ## corrected epsilon
                iH <- iOmega.adj.chol %*% solve(chol(crossprod(iOmega.adj) - iOmega.adj.chol %*% ls.Psi[[iC]] %*% iOmega.adj.chol)) %*% iOmega.adj.chol 
                epsilon.adj[[iC]] <- epsilon[[iC]] %*% iH

                ## derivative of the score regarding Y
                scoreY <- ls.dmu[[iC]] %*% iOmegaM1.adj
                for(iP in 1:n.varparam){ ## iP <- 1
                    scoreY[name.varparam[iP],] <- 2 * epsilon.adj[[iC]] %*% (iOmegaM1.adj %*% dOmega[[name.varparam[iP]]][iIndex.Omega,iIndex.Omega]  %*% iOmegaM1.adj)
                }
                ## leverage
                leverage[[iC]] <- colSums(iVcov.param %*% ls.dmu[[iC]] * scoreY) ## NOTE: dimensions of ls.dmu and scoreY matches even when there are missing values
                # same as
                # diag(t(ls.dmu[[iC]])  %*% iVcov.param %*% scoreY)
            }
            if(is.null(index.Omega)){ ## corrected sample size                
                n.corrected <- rep(n.cluster, n.endogenous) - Reduce("+",leverage)
            }
        }

      
        ## *** Step (i): expected information matrix
        iInfo <- .information2(dmu = dmu,
                               dOmega = dOmega,
                               Omega = Omega.adj,
                               OmegaM1 = OmegaM1.adj,
                               n.corrected = n.corrected,
                               leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                               grid.meanparam = grid.meanparam,
                               n.grid.meanparam = n.grid.meanparam,
                               grid.varparam = grid.varparam,
                               n.grid.varparam = n.grid.varparam,
                               name.param = name.param,
                               name.meanparam = name.meanparam,
                               name.varparam = name.varparam,
                               param2index = param2index, n.param = n.param)
        iVcov.param <- solve(iInfo)
        
        ## *** Update cv
        iIter <- iIter + 1
        iTol <- norm(Omega.adj-Omega_save, type = "F")
        ## cat("Omega.adj: ",Omega.adj," | n:",n.corrected," | iTol:",iTol,"\n")
    }

    ## ** Post processing
    dimnames(OmegaM1.adj) <- list(name.endogenous, name.endogenous)

    vcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if("try-error" %in% class(vcov.param)){
        errorMessage <- vcov.param
        vcov.param <- solve(iInfo)
        attr(vcov.param, "warning") <- errorMessage
    }
    dimnames(vcov.param) <- dimnames(iInfo)

    if(is.null(index.Omega)){
       leverage <- do.call(rbind,leverage)
       epsilon.adj <- do.call(rbind,epsilon.adj)        
    }
    
    ## ** Export    
    return(list(Omega = Omega.adj,
                OmegaM1 = OmegaM1.adj,
                vcov.param = vcov.param,
                leverage = leverage,
                n.corrected = n.corrected,
                epsilon = epsilon.adj,
                iter = iIter))
}


##----------------------------------------------------------------------
### adjustEstimate.R ends here
