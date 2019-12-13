### estimate2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: dec 13 2019 (17:26) 
##           By: Brice Ozenne
##     Update #: 898
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
#' @name estimate2
#' 
#' @keywords internal
.sscResiduals <- function(object, param, algorithm = "2"){

    algorithm <- match.arg(as.character(algorithm), choices = c("1","2"))

    ## ** current values
    param <- object$sCorrect$param
    Omega <- object$sCorrect$moment$Omega
    epsilon <- object$sCorrect$moment$residuals
    leverage <- object$sCorrect$leverage
    dmu <- object$sCorrect$dmoment$dmu
    dOmega <- object$sCorrect$dmoment$dOmega
    vcov.param <- object$sCorrect$vcov.param
    
    endogenous <- object$sCorrect$endogenous
    n.endogenous <- length(endogenous)
    n.cluster <- object$sCorrect$cluster$n.cluster
    param.mean <- object$sCorrect$skeleton$Uparam.mean
    param.var <- object$sCorrect$skeleton$Uparam.var
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
    
    for(iP in 1:n.pattern){ ## iP <- 1
        iY <- unique.pattern[iP,]
        
        for(iC in missing.pattern[[iP]]){ ## iC <- 1
            ## individual bias
            iPsi <- t(dmu[,iY,iC])  %*% vcov.muparam %*% dmu[,iY,iC]

            ## cumulated bias            
            Psi[iY,iY] <- Psi[iY,iY] + iPsi
            n.Psi[iY,iY] <- n.Psi[iY,iY] + 1
        }
    }

    ## take the average
    Psi[n.Psi>0] <- Psi[n.Psi>0]/n.Psi[n.Psi>0]
        
    ## ** Step (iii): compute corrected residuals and effective sample size
    if(algorithm == "2"){
        epsilon.adj <- .adjustResiduals(Omega = Omega,
                                        Psi = Psi,
                                        epsilon = epsilon,
                                        name.pattern = name.pattern,
                                        missing.pattern = missing.pattern,
                                        unique.pattern = unique.pattern,
                                        endogenous = endogenous,
                                        n.endogenous = n.endogenous,
                                        n.cluster = n.cluster)

        leverage.adj <- .leverage2(Omega = Omega,
                                   epsilon = epsilon.adj,
                                   dmu = dmu,
                                   dOmega = dOmega,
                                   name.pattern = name.pattern,
                                   missing.pattern = missing.pattern,
                                   unique.pattern = unique.pattern,
                                   endogenous = endogenous,
                                   n.endogenous = n.endogenous,
                                   param.var = param.var,
                                   n.param.var = length(param.var),
                                   n.cluster = n.cluster)

     
    }else{
        epsilon.adj <- epsilon
        leverage.adj <- leverage
    }
        
    ## ** Step (iv): bias-corrected residual covariance matrix
    Omega.adj <- Omega + Psi

    ## ** Step (v): bias-corrected variance parameters

    ## *** extract info
    n.endogenous <- NROW(Omega)
    df.param <- object$conditionalMoment$df.param
    
    index.matrix <- object$conditionalMoment$adjustMoment$index.matrix
    index.Psi <- object$conditionalMoment$adjustMoment$index.Psi
    A <- object$conditionalMoment$adjustMoment$A
    name.var <- object$conditionalMoment$adjustMoment$name.var
    n.rhs <- object$conditionalMoment$adjustMoment$n.rhs
    index.LambdaB <- object$conditionalMoment$adjustMoment$index.LambdaB
    name.endogenous <- object$conditionalMoment$adjustMoment$name.endogenous
    name.latent <- object$conditionalMoment$adjustMoment$name.latent
    
    skeleton <- object$conditionalMoment$skeleton

    param <- object$conditionalMoment$param
    Lambda <- object$conditionalMoment$value$Lambda
    iIB <- object$conditionalMoment$value$iIB
    iIB.Lambda <- object$conditionalMoment$value$iIB.Lambda
    dLambda <- object$conditionalMoment$dMoment.init$dLambda
    dB <- object$conditionalMoment$dMoment.init$dB

    ## *** right hand side of the equation
    eq.rhs <- Omega[index.matrix$index]

    ## *** left hand side of the equation
    if(NROW(index.Psi)>0){
        n.index.Psi <- NROW(index.Psi)
        n.latent <- NROW(skeleton$Psi)        
        Z <- iIB %*% Lambda

        ## A = t(Z) Psi Z + Sigma
        ## (t(Z) Psi Z)_{ij} = \sum_{k,l} Z_{k,i} Psi_{k,l} Z_{l,j}
        for(iIndex in 1:n.rhs){ # iIndex <- 1
            iRow <- index.matrix[iIndex,"row"]
            iCol <- index.matrix[iIndex,"col"]
            for(iPsi in 1:n.index.Psi){
                iRowPsi <- index.Psi[iPsi,"row"]
                iColPsi <- index.Psi[iPsi,"col"]
                A[iIndex,skeleton$Psi[iRowPsi,iColPsi]] <- A[iIndex,skeleton$Psi[iRowPsi,iColPsi]] + Z[iRowPsi,iRow]*Z[iColPsi,iCol]
            }
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
    object$conditionalMoment$param[name.var] <- setNames(iSolution, name.var)

    ## ** Step (vi-vii): update derivatives and information matrix (performed by .init_sCorrect) in the level above
        
    ## ** Export
    return(newparam)
}



##----------------------------------------------------------------------
### estimate2.R ends here
