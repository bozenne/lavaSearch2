### sCorrect-sscResiduals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: jan  8 2020 (11:15) 
##           By: Brice Ozenne
##     Update #: 989
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .init_sscResiduals
.init_sscResiduals <- function(object){

    out <- list()

    ## ** extract info
    endogenous <- object$sCorrect$endogenous
    latent <- object$sCorrect$latent
    
    type <- object$sCorrect$skeleton$type
    type <- type[!is.na(type$param),]
    type <- type[type$detail %in% c("Sigma_var","Sigma_cov","Psi_var","Psi_cov","sigma2","sigma2k","cor"),]

    Omega <- object$sCorrect$moment$Omega
    name.var <- unique(type$param)
    n.var <- length(name.var)

    ## ** subset residual variance-covariance
    index.upper.tri <- data.frame(index = which(upper.tri(Omega, diag = TRUE)),
                                  which(upper.tri(Omega, diag = TRUE), arr.ind = TRUE)
                                  )
    name.rhs <- paste(endogenous[index.upper.tri[,"row"]],
                      lava.options()$symbols[2],
                      endogenous[index.upper.tri[,"col"]],
                      sep = "")
    n.rhs <- length(name.rhs)

    ## ** Design matrix
    A <- matrix(0, nrow = n.rhs, ncol = n.var,
                    dimnames = list(name.rhs, name.var))

    ## *** Sigma_var and Sigma_cov
    if(any(type$detail %in% c("Sigma_var","Sigma_cov"))){
        type.Sigma <- type[type$detail %in% c("Sigma_var","Sigma_cov"),,drop=FALSE]
        for(iRow in 1:NROW(type.Sigma)){ ## iRow <- 1 
            A[paste0(type.Sigma[iRow,"Y"],"~~",type.Sigma[iRow,"X"]),type.Sigma[iRow,"param"]] <- 1
        }
    }
    ## *** sigma2, sigma2k, and cor
    if(any(type$detail %in% c("sigma2","sigma2k","cor"))){
        stop("Not yet implemented")
    }
    
    ## *** Psi_var and Psi_cov
    if(any(type$detail %in% c("Psi_var","Psi_cov"))){
        type.Psi <- type[type$detail %in% c("Psi_var","Psi_cov"),,drop=FALSE]
        index.Psi <- cbind(row = match(type.Psi$X, latent),
                           col = match(type.Psi$Y, latent))
        rownames(index.Psi) <- type.Psi$param
    }else{
        index.Psi <- NULL
    }

    return(list(type = "residuals",
                param0 = object$sCorrect$param,
                Omega0 = object$sCorrect$moment$Omega,
                residuals0 = object$sCorrect$residuals,
                index.upper.tri = index.upper.tri,
                name.rhs = name.rhs,
                name.var = name.var,
                A = A,
                index.Psi = index.Psi
                ))
}

## * .sscResiduals
#' @title Compute Bias Corrected Quantities.
#' @description Compute bias corrected residuals variance covariance matrix
#' and information matrix.
#' Also provides the leverage values and corrected sample size when adjust.n is set to TRUE.
#' @name estimate2
#' 
#' @keywords internal
.sscResiduals <- function(object, param, algorithm = "2"){
    algorithm <- match.arg(as.character(algorithm), choices = c("1","2"))

    ## ** current values
    Omega0 <- object$sCorrect$ssc$Omega0 ## non bias corrected value of Omega
    param0 <- object$sCorrect$ssc$param0 ## non bias corrected value of the model parameters
    residuals0 <- object$sCorrect$ssc$residuals0 ## non bias corrected value of the model residuals
    
    epsilon <- object$sCorrect$residuals
    leverage <- object$sCorrect$leverage
    dmu <- object$sCorrect$dmoment$dmu
    dOmega <- object$sCorrect$dmoment$dOmega
    vcov.param <- object$sCorrect$vcov.param
    
    endogenous <- object$sCorrect$endogenous
    n.endogenous <- length(endogenous)
    n.cluster <- object$sCorrect$cluster$n.cluster
    param.mean <- object$sCorrect$skeleton$Uparam.mean
    param.var <- object$sCorrect$skeleton$Uparam.var
    param.hybrid <- intersect(param.mean,param.var)
    missing.pattern <- object$sCorrect$missing$pattern
    name.pattern <- object$sCorrect$missing$name.pattern
    n.pattern <- length(name.pattern)
    unique.pattern <- object$sCorrect$missing$unique.pattern

    ## ** Step (i-ii) compute individual and average bias
    dmu <- aperm(abind::abind(dmu, along = 3), perm = c(3,2,1))
    vcov.muparam <- vcov.param[param.mean,param.mean,drop=FALSE]

    Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                  dimnames = list(endogenous, endogenous))
    n.Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                    dimnames = list(endogenous, endogenous))

    ls.Psi <- vector(mode = "list", length = n.cluster)
    for(iP in 1:n.pattern){ ## iP <- 1
        iY <- unique.pattern[iP,]
        
        for(iC in missing.pattern[[iP]]){ ## iC <- 1
            ## individual bias
            iPsi <- t(dmu[,iY,iC])  %*% vcov.muparam %*% dmu[,iY,iC]
            ls.Psi[[iC]] <- iPsi
            ## cumulated bias            
            Psi[iY,iY] <- Psi[iY,iY] + iPsi
            n.Psi[iY,iY] <- n.Psi[iY,iY] + 1
        }
    }

    ## take the average
    Psi[n.Psi>0] <- Psi[n.Psi>0]/n.Psi[n.Psi>0]
    
    ## ** Step (iii): compute corrected residuals and effective sample size
    if(algorithm == "2"){
        epsilon.adj <- .adjustResiduals(Omega = Omega0,
                                        Psi = Psi,
                                        epsilon = residuals0,
                                        name.pattern = name.pattern,
                                        missing.pattern = missing.pattern,
                                        unique.pattern = unique.pattern,
                                        endogenous = endogenous,
                                        n.cluster = n.cluster)

        leverage.adj <- .leverage2(Omega = Omega0,
                                   epsilon = epsilon.adj,
                                   dmu = dmu,
                                   dOmega = dOmega,
                                   vcov.param = vcov.param,
                                   name.pattern = name.pattern,
                                   missing.pattern = missing.pattern,
                                   unique.pattern = unique.pattern,
                                   endogenous = endogenous,
                                   n.endogenous = n.endogenous,
                                   param = object$sCorrect$skeleton$Uparam,
                                   param.mean = param.mean,
                                   param.hybrid = param.hybrid,
                                   n.cluster = n.cluster)

     
    }else{
        epsilon.adj <- epsilon
        leverage.adj <- leverage
    }
        
    ## ** Step (iv): bias-corrected residual covariance matrix
    Omega.adj <- Omega0 + Psi

    ## ** Step (v): bias-corrected variance parameters

    ## *** right hand side of the equation
    index.upper.tri <- object$sCorrect$ssc$index.upper.tri[,"index"]
    eq.rhs <- setNames(Omega.adj[index.upper.tri],
                       object$sCorrect$ssc$name.rhs)

    ## *** left hand side of the equation
    A <- object$sCorrect$ssc$A
    index.Psi <- object$sCorrect$ssc$index.Psi

    if(NROW(index.Psi)>0){
        Z <- object$sCorrect$moment$iIB %*% object$sCorrect$moment$Lambda
        tZ <- t(Z)
        n.index.Psi <- NROW(index.Psi)
    
        ## A = t(Z) Psi Z + Sigma
        ## (t(Z) Psi Z)_{ij} = \sum_{k,l} Z_{k,i} Psi_{k,l} Z_{l,j}
        ## (t(Z) Psi Z)_{ij} regarding param_(k,l) = Z_{k,i} Z_{l,j}
        for(iPsi in 1:n.index.Psi){ ## iPsi <- 1
            iNamePsi <- rownames(index.Psi)[iPsi]
            iRowPsi <- index.Psi[iPsi,"row"]
            iColPsi <- index.Psi[iPsi,"col"]
            A[,iNamePsi] <- A[,iNamePsi] + (tZ[,index.Psi[iPsi,"row"]] %o% Z[index.Psi[iPsi,"col"],])[index.upper.tri]
        }
        
    }

    ## *** solve equation
    ## microbenchmark::microbenchmark(svd = {asvd <- svd(A) ; asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs;},
    ## qr = qr.coef(qr(A), eq.rhs),
    ## Rcpp = OLS_cpp(A, eq.rhs),
    ## RcppTry = try(OLS_cpp(A, eq.rhs)[,1], silent = TRUE),
    ## Rcpp2 = OLS2_cpp(A, eq.rhs),
    ## OLS1 = solve(crossprod(A), crossprod(A, eq.rhs)),
    ## OLS2 = solve(t(A) %*% A) %*% t(A) %*% eq.rhs,
    ## OLS_stats = stats::lsfit(x = A, y = eq.rhs),
    ## OLS_LINPACK = .Call(stats:::C_Cdqrls, x = A, y = eq.rhs, tolerance = 1e-7, FALSE)$coefficients, times = 500)
    if(lava.options()$method.estimate2=="svd"){
        asvd <- svd(A)
        iSolution <- try((asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs)[,1], silent = TRUE)
    }else if(lava.options()$method.estimate2=="ols"){
        iSolution <- try(OLS_cpp(A, eq.rhs)[,1], silent = TRUE)
    }else{
        stop("unknown OLS methods \n")
    }
    
    if(inherits(iSolution, "try-error")){
        if(abs(det(t(A) %*% A)) <  1e-10){            
            stop("Singular matrix: cannot update the estimates \n")
        }else{
            stop(iSolution)
        }
    }

    ## ** update parameters in conditional moments
    ## iSolution - param[object$sCorrect$ssc$name.var]
    param0[object$sCorrect$ssc$name.var] <- iSolution

    ## ** Step (vi-vii): update derivatives and information matrix (performed by .init_sCorrect) in the parent function
        
    ## ** Export
    attr(param0,"leverage") <- leverage.adj
    attr(param0,"residuals") <- epsilon.adj
    attr(param0,"Psi") <- Psi
    return(param0)
}



##----------------------------------------------------------------------
### sCorrect-sscResiduals.R ends here
