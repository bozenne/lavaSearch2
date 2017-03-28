### estimate2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: mar 28 2017 (17:25) 
##           By: Brice Ozenne
##     Update #: 651
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * estimate2
#' @title Compute Bias Corrected Quantities.
#' @description Compute bias corrected residuals variance covariance matrix
#' and information matrix.
#' Also provides the leverage values and corrected sample size when adjust.n is set to TRUE.
#' 
#' @keywords internal
.estimate2 <- function(object, epsilon, n.cluster,
                       name.param, name.endogenous, name.meanparam, name.varparam,
                       index.Omega,
                       adjust.Omega, adjust.n, tol, n.iter, trace){

    ## ** Prepare
    Omega <- object$conditionalMoment$Omega
    dmu <- object$conditionalMoment$dmu
    dOmega <- object$conditionalMoment$dOmega
    
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

    ## check low diagonal
    name2num <- setNames(1:n.param,name.param)
    if(!all(name2num[grid.meanparam[,1]]-name2num[grid.meanparam[,2]]>=0)){
        stop("Incorrect allocation of the computation of the information matrix (mean parameter) \n")
    }
    name2num <- setNames(1:n.param,name.param)
    if(!all(name2num[grid.varparam[,1]]-name2num[grid.varparam[,2]]>=0)){
        stop("Incorrect allocation of the computation of the information matrix (variance parameter) \n")
    }
    ##
    
    param2index <- setNames(1:n.param, name.param)

    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, name.endogenous))
    ls.dmu <- vector(mode = "list", length = n.cluster)
    for(iC in 1:n.cluster){ # iC <- 1
        if(is.null(index.Omega)){            
            leverage[iC,] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = n.endogenous,
                                   dimnames = list(name.param, name.endogenous))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu[name.meanparam],function(x){x[iC,]}))
        }else{
            leverage[iC,index.Omega[[iC]]] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = length(index.Omega[[iC]]),
                                   dimnames = list(name.param, name.endogenous[index.Omega[[iC]]]))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu[name.meanparam],function(x){x[iC,index.Omega[[iC]]]}))
        }        
    }
    
    ## ** Initialisation (i.e. first iteration without correction)
    if(any(eigen(Omega)$value<=0)){
        stop("the residual variance-covariance matrix is not positive definite \n")
    }

    if(is.null(index.Omega)){
        n.corrected <- rep(n.cluster, n.endogenous)
    }else{
        n.corrected <- NULL
    }
    ls.Psi <- vector(mode = "list", length = n.cluster)

    Omega.adj <- Omega
    if(!adjust.n){
       epsilon.adj <- epsilon
    }

    if(trace>0){
        cat("* Reconstruct estimated information matrix ")
    }

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
                           name.meanparam = name.meanparam,
                           name.varparam = name.varparam,
                           param2index = param2index, n.param = n.param)
    iVcov.param <- chol2inv(chol(iInfo))
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
        Omega_save <- Omega
        iOmega.adj <- Omega.adj
    }else{
        iIter <- Inf
        iTol <- -Inf        
    }
    
    while(iIter < n.iter & iTol > tol){
        if(trace>0){
            cat("*")
        }

        ## *** Step (i-ii): compute individual bias, expected bias
        Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                      dimnames = list(name.endogenous, name.endogenous))
        M.countCluster <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                                 dimnames = list(name.endogenous, name.endogenous))
        for(iC in 1:n.cluster){
            ## individual bias
            ls.Psi[[iC]] <- t(ls.dmu[[iC]])  %*% iVcov.param %*% ls.dmu[[iC]]
            ## cumulated bias            
            if(is.null(index.Omega)){
                Psi <- Psi + ls.Psi[[iC]]
                M.countCluster <- M.countCluster + 1
            }else{
                Psi[index.Omega[[iC]],index.Omega[[iC]]] <- Psi[index.Omega[[iC]],index.Omega[[iC]]] + ls.Psi[[iC]]
                M.countCluster[index.Omega[[iC]],index.Omega[[iC]]] <- M.countCluster[index.Omega[[iC]],index.Omega[[iC]]] + 1
            }
        }

        ## update
        for(iPsi in 1:length(Psi)){
            if(M.countCluster[iPsi]>0){
                Psi[iPsi] <- Psi[iPsi]/M.countCluster[iPsi]
            }
        }
        
        ## *** Step (iii): compute leverage
        if(adjust.n){
            res.tempo <- .calcLeverage(Omega = Omega.adj,
                                       ls.Psi = ls.Psi,
                                       epsilon = epsilon,
                                       ls.dmu = ls.dmu,
                                       dOmega = dOmega,
                                       vcov.param = iVcov.param,
                                       index.Omega = index.Omega,
                                       name.endogenous = name.endogenous,
                                       n.endogenous = n.endogenous,
                                       name.varparam = name.varparam,
                                       n.varparam = n.varparam,
                                       n.cluster = n.cluster)
            leverage <- res.tempo$leverage
            epsilon.adj <- res.tempo$residuals
            n.corrected <- rep(n.cluster, n.endogenous) - colSums(leverage, na.rm = TRUE)
        }
        
        ## *** Step (iv): correct residual covariance matrix and estimates
        if(adjust.Omega){
            ## corrected residual covariance variance
            Omega.adj <- Omega + Psi
            
            ## correct estimates
            object$conditionalMoment <- .adjustMoment(object, Omega = Omega.adj)
            dOmega <- object$conditionalMoment$dOmega
            ## conditionalMoment.adj$param - coef(object)
           
        }

        ## *** Step (v): expected information matrix
        iInfo <- .information2(dmu = dmu,
                               dOmega = dOmega,
                               Omega = Omega.adj,
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
        iVcov.param <- chol2inv(chol(iInfo))
        
        ## *** Update cv
        iIter <- iIter + 1
        iTol <- norm(Omega.adj-Omega_save, type = "F")
        Omega_save <- Omega.adj
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

    vcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if("try-error" %in% class(vcov.param)){
        errorMessage <- vcov.param
        vcov.param <- solve(iInfo)
        attr(vcov.param, "warning") <- errorMessage
    }
    dimnames(vcov.param) <- dimnames(iInfo)

    ## update object
    object$conditionalMoment$Omega <- Omega.adj
    object$dVcov <- list(param = object$conditionalMoment$param,
                         score = NULL,
                         vcov.param = vcov.param,
                         dVcov.param = NULL,
                         Omega = Omega.adj,
                         residuals = epsilon.adj,
                         leverage = leverage,
                         n.corrected = rep(n.cluster, n.endogenous) - colSums(leverage, na.rm = TRUE),
                         opt = list(objective = iTol, iterations = iIter, convergence = (iTol <= tol)))

    ## ** Export
    return(object)
}

## * .calcLeverage
.calcLeverage <- function(Omega, ls.Psi, epsilon, ls.dmu, dOmega, vcov.param,
                          index.Omega,
                          name.endogenous, n.endogenous, name.varparam, n.varparam, n.cluster){

    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, name.endogenous))
    epsilon.adj <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                          dimnames = list(NULL, name.endogenous))
    
    ## ** symmetric square root.
    ## Omega.adj.chol %*% Omega.adj.chol - Omega.adj
    ## Omega.adj.chol %*% Omega.adj %*% Omega.adj.chol - Omega.adj %*% Omega.adj
    if(is.null(index.Omega)){
        iOmega <- Omega
        iOmega.chol <- matrixPower(iOmega, symmetric = TRUE, power = 1/2)
        iOmegaM1 <- chol2inv(iOmega.chol)
        iIndex.Omega <- 1:n.endogenous
    }
            
    for(iC in 1:n.cluster){                 # iC <- 1
        if(!is.null(index.Omega)){
            iIndex.Omega <- index.Omega[[iC]]
            iOmega <- Omega[iIndex.Omega,iIndex.Omega,drop=FALSE]
            iOmega.chol <- matrixPower(iOmega, symmetric = TRUE, power = 1/2)
            iOmegaM1 <- chol2inv(iOmega.chol)
        }
        ## corrected epsilon
        H <- iOmega %*% iOmega - iOmega.chol %*% ls.Psi[[iC]] %*% iOmega.chol
        ## iH <- matrixPower(H, symmetric = TRUE, power = -1/2)
        iH <- tryCatch(matrixPower(H, symmetric = TRUE, power = -1/2), warning = function(w){w})
        if(inherits(iH,"warning")){
            stop("Cannot compute the adjusted residuals \n",
                 "Estimated bias too large compared to the estimated variance-covariance matrix \n",
                 "Consider setting argument \'adjust.n\' to FALSE when calling sCorrect \n")
        }
        ## Omega.adj.chol %*% iH %*% Omega.adj.chol %*% (iOmega.adj - ls.Psi[[iC]]) %*% Omega.adj.chol %*% iH %*% Omega.adj.chol - Omega.adj
        epsilon.adj[iC,iIndex.Omega] <- epsilon[iC,iIndex.Omega] %*% iOmega.chol %*% iH %*% iOmega.chol

        ## derivative of the score regarding Y
        scoreY <- ls.dmu[[iC]] %*% iOmegaM1
        for(iP in 1:n.varparam){ ## iP <- 1
            scoreY[name.varparam[iP],] <- scoreY[name.varparam[iP],] + 2 * epsilon.adj[iC,iIndex.Omega] %*% (iOmegaM1 %*% dOmega[[name.varparam[iP]]][iIndex.Omega,iIndex.Omega]  %*% iOmegaM1)
        }
        ## leverage
        leverage[iC,iIndex.Omega] <- colSums(vcov.param %*% ls.dmu[[iC]] * scoreY) ## NOTE: dimensions of ls.dmu and scoreY matches even when there are missing values
                                        # same as
                                        # diag(t(ls.dmu[[iC]])  %*% iVcov.param %*% scoreY)
    }

    return(list(leverage = leverage,
                residuals = epsilon.adj))            
}

## * .adjustMoment
`.adjustMoment` <-
    function(object, ...) UseMethod(".adjustMoment")

## * .adjustMoment.lvmfit

.adjustMoment.lvmfit <- function(object, Omega){

    ## ** extract info
    n.endogenous <- NROW(Omega)
    df.param <- object$conditionalMoment$df.param
    
    index.matrix <- object$conditionalMoment$adjustMoment$index.matrix
    index.Psi <- object$conditionalMoment$adjustMoment$index.Psi
    A <- object$conditionalMoment$adjustMoment$A
    name.var <- object$conditionalMoment$adjustMoment$name.var
    n.rhs <- object$conditionalMoment$adjustMoment$n.rhs

    Psi <- object$conditionalMoment$skeleton$Psi
    Lambda <- object$conditionalMoment$value$Lambda
    iIB <- object$conditionalMoment$value$iIB
    iIB.Lambda <- object$conditionalMoment$value$iIB.Lambda
    dLambda <- object$conditionalMoment$dMoment.init$dLambda
    dB <- object$conditionalMoment$dMoment.init$dB
    
    ## ** right hand side of the equation
    eq.rhs <- Omega[index.matrix$index]
    
    ## ** left hand side of the equation
    if(NROW(index.Psi)>0){
        n.index.Psi <- NROW(index.Psi)
        n.latent <- NROW(Psi)        
        Z <- iIB %*% Lambda

        ## A = t(Z) Psi Z + Sigma
        ## (t(Z) Psi Z)_{ij} = \sum_{k,l} Z_{k,i} Psi_{k,l} Z_{l,j}
        for(iIndex in 1:n.rhs){ # iIndex <- 1
            iRow <- index.matrix[iIndex,"row"]
            iCol <- index.matrix[iIndex,"col"]
            for(iPsi in 1:n.index.Psi){
                iRowPsi <- index.Psi[iPsi,"row"]
                iColPsi <- index.Psi[iPsi,"col"]
                A[iIndex,Psi[iRowPsi,iColPsi]] <- A[iIndex,Psi[iRowPsi,iColPsi]] + Z[iRowPsi,iRow]*Z[iColPsi,iCol]
            }
        }
    }
    
    ## ** solve equation
    asvd <- svd(A)
    if(any(abs(asvd$d) < .Machine$double.eps ^ 0.5)){
        stop("Singular matrix: cannot update the estimates \n")
    }
    object$conditionalMoment$param[name.var] <- setNames(as.double(asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs),
                                                         name.var)
    

    ## ** update conditional moments
    ## *** Sigma
    if(length(object$conditionalMoment$skeleton$Sigma)>0){
        index.update <- which(!is.na(object$conditionalMoment$skeleton$Sigma))
        object$conditionalMoment$value$Sigma[index.update] <- object$conditionalMoment$param[object$conditionalMoment$skeleton$Sigma[index.update]]
    }

    ## *** Psi
    if(length(object$conditionalMoment$skeleton$Psi)>0){
        index.update <- which(!is.na(object$conditionalMoment$skeleton$Psi))
        object$conditionalMoment$value$Psi[index.update] <- object$conditionalMoment$param[object$conditionalMoment$skeleton$Psi[index.update]]

         object$conditionalMoment$value$tLambda.tiIB.Psi.iIB <- t(Lambda) %*% t(iIB) %*% object$conditional$value$Psi %*% iIB    
    }

    ## ** update first derivative of the conditional variance
    ## only Lambda and B
    name.vcovparam <- names(object$conditionalMoment$dOmega)
    type.vcovparam <- df.param[match(name.vcovparam,df.param$name),"detail"]

    indexUpdate <- which(type.vcovparam %in% c("Lambda","B"))
    name2update <- name.vcovparam[indexUpdate]
    type2update <- type.vcovparam[indexUpdate]
    n.update <- length(type2update)

    if(n.update>0){

        for(iP in 1:n.update){ # iP <- 1
            iType <- type2update[iP]
            iName <- name2update[iP]

            if(iType == "Lambda"){
                object$conditionalMoment$dOmega[[iName]] <- object$conditionalMoment$value$tLambda.tiIB.Psi.iIB %*% dLambda[[iName]]
                object$conditionalMoment$dOmega[[iName]] <- object$conditionalMoment$dOmega[[iName]] + t(object$conditionalMoment$dOmega[[iName]])
            }else if(iType == "B"){
                object$conditionalMoment$dOmega[[iName]] <- object$conditionalMoment$value$tLambda.tiIB.Psi.iIB %*% dB[[iName]] %*% iIB.Lambda
                object$conditionalMoment$dOmega[[iName]] <- object$conditionalMoment$dOmega[[iName]] + t(object$conditionalMoment$dOmega[[iName]])
            }

        }
    }
    
    ## ** export
    return(object$conditionalMoment)
}

##----------------------------------------------------------------------
### estimate2.R ends here
