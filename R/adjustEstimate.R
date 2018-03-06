### adjustEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: mar  6 2018 (11:44) 
##           By: Brice Ozenne
##     Update #: 366
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
                           index.Omega,
                           adjust.Omega, adjust.n, tol, n.iter, trace){

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
        n.corrected <- rep(n.cluster, n.endogenous)        
    }else{
        n.corrected <- NULL        
    }

    ls.Psi <- vector(mode = "list", length = n.cluster)
    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, name.endogenous))
    ls.dmu <- vector(mode = "list", length = n.cluster)
    for(iC in 1:n.cluster){ # iC <- 1
        if(is.null(index.Omega)){
            leverage[iC,] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = n.endogenous,
                                   dimnames = list(name.param, name.endogenous))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu,function(x){x[iC,]}))
        }else{
            leverage[iC,index.Omega[[iC]]] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = length(index.Omega[[iC]]),
                                   dimnames = list(name.param, name.endogenous[index.Omega[[iC]]]))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu,function(x){x[iC,index.Omega[[iC]]]}))
        }        
    }
    ## ** Initialisation (i.e. first iteration without correction)
    Omega.adj <- Omega
    OmegaM1.adj <- chol2inv(chol(Omega.adj))

    if(trace>0){
        cat("* Reconstruct estimated information matrix ")
    }
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
    if(trace>0){
        cat("- done \n")
    }
    
    ## ** Loop    
    if(adjust.Omega || adjust.n){
        if(trace>0){
            cat("* iterative small sample correction: ")
        }
        iIter <- 0
        iTol <- Inf
    }else{
        iIter <- Inf
        iTol <- -Inf
    }

    while(iIter < n.iter & iTol > tol){
        if(trace>0){
            cat("*")
        }
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
            ## symmetric square root.
            Omega.adj.chol <- matrixPower(Omega.adj, symmetric = TRUE, power = 1/2)
            ## Omega.adj.chol %*% Omega.adj.chol - Omega.adj
            ## Omega.adj.chol %*% Omega.adj %*% Omega.adj.chol - Omega.adj %*% Omega.adj
            if(is.null(index.Omega)){
                iOmega.adj <- Omega.adj
                iOmegaM1.adj <- OmegaM1.adj
                iOmega.adj.chol <- Omega.adj.chol
                iIndex.Omega <- 1:n.endogenous
            }
            
            for(iC in 1:n.cluster){                 # iC <- 1
                if(!is.null(index.Omega)){
                    iOmega.adj <- Omega.adj[index.Omega[[iC]],index.Omega[[iC]],drop=FALSE]
                    iOmegaM1.adj <- OmegaM1.adj[index.Omega[[iC]],index.Omega[[iC]],drop=FALSE]
                    iOmega.adj.chol <- Omega.adj.chol[index.Omega[[iC]],index.Omega[[iC]],drop=FALSE]
                    iIndex.Omega <- index.Omega[[iC]]
                }
                ## corrected epsilon
                iH <- matrixPower(iOmega.adj %*% iOmega.adj - iOmega.adj.chol %*% ls.Psi[[iC]] %*% iOmega.adj.chol,
                                  symmetric = TRUE, power = -1/2)
                ## Omega.adj.chol %*% iH %*% Omega.adj.chol %*% (iOmega.adj - ls.Psi[[iC]]) %*% Omega.adj.chol %*% iH %*% Omega.adj.chol - Omega.adj
                epsilon.adj[iC,iIndex.Omega] <- epsilon[iC,iIndex.Omega] %*% iOmega.adj.chol %*% iH %*% iOmega.adj.chol
                ## derivative of the score regarding Y
                scoreY <- ls.dmu[[iC]] %*% iOmegaM1.adj
                for(iP in 1:n.varparam){ ## iP <- 1
                    scoreY[name.varparam[iP],] <- scoreY[name.varparam[iP],] + 2 * epsilon.adj[iC,iIndex.Omega] %*% (iOmegaM1.adj %*% dOmega[[name.varparam[iP]]][iIndex.Omega,iIndex.Omega]  %*% iOmegaM1.adj)
                }
                ## leverage
                leverage[iC,iIndex.Omega] <- colSums(iVcov.param %*% ls.dmu[[iC]] * scoreY) ## NOTE: dimensions of ls.dmu and scoreY matches even when there are missing values
                # same as
                # diag(t(ls.dmu[[iC]])  %*% iVcov.param %*% scoreY)
            }
            if(is.null(index.Omega)){ ## corrected sample size                
                n.corrected <- rep(n.cluster, n.endogenous) - colSums(leverage)
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
    if(!is.infinite(iIter)){

        if(iTol > tol){
            warning("small sample correction did not reach convergence after ",iIter," iterations \n")

            if(trace>0){
                cat(" - incomplete \n")
            }
        }else{
            if(trace>0){
                cat(" - done \n")
            }
        }
        
    }

    dimnames(OmegaM1.adj) <- list(name.endogenous, name.endogenous)

    vcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if("try-error" %in% class(vcov.param)){
        errorMessage <- vcov.param
        vcov.param <- solve(iInfo)
        attr(vcov.param, "warning") <- errorMessage
    }
    dimnames(vcov.param) <- dimnames(iInfo)

    if(!is.null(index.Omega)){
        n.corrected <- colSums(leverage, na.rm = TRUE)
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
