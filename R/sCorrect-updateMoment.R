### sCorrect-updateMoment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 10 2019 (09:58) 
## Version: 
## Last-Updated: jan 15 2020 (15:33) 
##           By: Brice Ozenne
##     Update #: 189
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * updateMoment
#' @rdname updateMoment
updateMoment <- function(skeleton, value, toUpdate,
                         name.pattern, unique.pattern,
                         param, endogenous, latent, n.cluster){
    if(lava.options()$debug){cat("updateMoment \n")}

    n.endogenous <- length(endogenous)
    n.latent <- length(latent)
    
    ## ** Update with the current values
    name.update  <- names(toUpdate[toUpdate==TRUE])
    if(length(name.update)>0){
        for(iUpdate in name.update){ ## iUpdate <- "Sigma"
            if(iUpdate == "SigmaValue"){
                index.update <- which(!is.na(skeleton$SigmaParam))
                skeleton$SigmaValue[index.update] <- param[skeleton$SigmaParam[index.update]]
                value$Sigma <- apply(skeleton$SigmaValue, MARGIN = 1:2, FUN = prod)
                attr(value$Sigma,"detail") <- skeleton$SigmaValue
            }else{
                index.update <- which(!is.na(skeleton[[iUpdate]]))
                value[[iUpdate]][index.update] <- param[skeleton[[iUpdate]][index.update]]
                skeleton$Sigma
            }
        }
    }

    ## ** Pre-compute relevant quantities
    if(n.latent>0){
        ## alpha + X\Gamma
        value$alpha.XGamma <- matrix(0, nrow = n.cluster, ncol = n.latent,
                                     dimnames = list(NULL, latent))
        if("alpha" %in% names(value)){
            value$alpha.XGamma <- value$alpha.XGamma + skeleton$Xalpha %o% value$alpha
        }
        if("Gamma" %in% names(value)){
            value$alpha.XGamma <- value$alpha.XGamma + do.call(cbind,lapply(latent, function(iL){skeleton$XGamma[[iL]] %*% value$Gamma[,iL]}))
        }

        ## (I-B)^{-1}
        if("B" %in% names(value)){
            value$iIB <- solve(diag(1, nrow = n.latent, ncol = n.latent) - value$B)
        }else{
            value$iIB <- diag(1, nrow = n.latent, ncol = n.latent)
            dimnames(value$iIB) <- list(latent,latent)
        }

        ## (alpha + X\Gamma) (I-B)^{-1}
        value$alpha.XGamma.iIB <- value$alpha.XGamma %*% value$iIB
        
        ## (I-B)^{-1} \Lambda
        value$iIB.Lambda <-  value$iIB %*% value$Lambda
        value$tLambda.tiIB <-  t(value$iIB.Lambda)

        ## \Psi (I-B)^{-1}
        value$Psi.iIB <- value$Psi %*% value$iIB

        ## (I-B)^{-t} \Psi (I-B)^{-1}
        value$tiIB.Psi.iIB <-  t(value$iIB) %*% value$Psi

        ## \Lambda^t (I-B)^{-t} \Psi (I-B)^{-1}
        value$tLambda.tiIB.Psi.iIB <- t(value$iIB.Lambda) %*% value$Psi.iIB

        ## (I-B)^{-t} \Psi (I-B)^{-1} \Lambda
        value$tiIB.Psi.iIB.Lambda <- t(value$tLambda.tiIB.Psi.iIB) 
    }

    ## ** Compute mean
    value$mu <- matrix(0, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL,endogenous))

    if("nu" %in% names(value)){
        value$mu <- value$mu + sweep(skeleton$Xnu, MARGIN = 2, FUN = "*", STATS = value$nu)
    }
    if("K" %in% names(value)){
        value$mu <- value$mu + do.call(cbind,lapply(endogenous, function(iE){skeleton$XK[[iE]] %*% value$K[,iE]})) ## iE <- endogenous[1]
    }
    if(n.latent>0){
        value$mu <- value$mu + value$alpha.XGamma %*% value$iIB.Lambda
    }

    ## ** Compute variance
    value$Omega <- matrix(0, nrow = n.endogenous, ncol = n.endogenous, 
                          dimnames = list(endogenous,endogenous))

    if("Sigma" %in% names(value)){
        value$Omega <- value$Omega + value$Sigma
    }
    if("Psi" %in% names(value)){
        value$Omega <- value$Omega + value$tLambda.tiIB.Psi.iIB %*% value$Lambda
    }

    value$Omega.missing.pattern <- lapply(1:length(name.pattern), function(iM){ ## iM <- 1
        iIndex <- which(unique.pattern[iM,]==1)
        return(value$Omega[iIndex,iIndex,drop=FALSE])
    })
    names(value$Omega.missing.pattern) <- name.pattern
    value$OmegaM1.missing.pattern <- lapply(value$Omega.missing.pattern, solve)

    ## ** Export
    return(value)
}

## * updateDMoment
#' @rdname updateMoment
updateDMoment <- function(moment, skeleton, param){
    if(lava.options()$debug){cat("updateDMoment \n")}

    ## ** import information
    dmu <- skeleton$dmu.dparam
    dOmega <- skeleton$dOmega.dparam

    iIB.Lambda <- moment$iIB.Lambda
    tLambda.tiIB <- moment$tLambda.tiIB
    alpha.XGamma.iIB <- moment$alpha.XGamma.iIB
    tiIB.Psi.iIB.Lambda <- moment$tiIB.Psi.iIB.Lambda
    tLambda.tiIB.Psi.iIB <- moment$tLambda.tiIB.Psi.iIB
    Sigma <- moment$Sigma
    attr(Sigma,"detail") <- NULL
        
    ## ** Compute partial derivative regarding the mean
    if("alpha" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$alpha)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dmu[[iParam]] <- skeleton$dmat.dparam$alpha[[iParam]] %*% iIB.Lambda
        }
    }

    if("Gamma" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$Gamma)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dmu[[iParam]] <- skeleton$dmat.dparam$Gamma[[iParam]] %*% iIB.Lambda
        }
    }

    if("Lambda" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$Lambda)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dmu[[iParam]] <- alpha.XGamma.iIB %*% skeleton$dmat.dparam$Lambda[[iParam]]
        }
    }

    if("B" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$B)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dmu[[iParam]] <- alpha.XGamma.iIB %*% skeleton$dmat.dparam$B[[iParam]] %*% iIB.Lambda
        }
    }
    
    ## ** Compute partial derivative regarding the variance
    if("Lambda" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$Lambda)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dOmega[[iParam]] <- t(skeleton$dmat.dparam$Lambda[[iParam]]) %*% tiIB.Psi.iIB.Lambda + tLambda.tiIB.Psi.iIB %*% skeleton$dmat.dparam$Lambda[[iParam]]
        }
    }

    if("B" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$B)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dOmega[[iParam]] <- t(iIB.Lambda) %*% t(skeleton$dmat.dparam$B[[iParam]]) %*% tiIB.Psi.iIB.Lambda + tLambda.tiIB.Psi.iIB %*% skeleton$dmat.dparam$B[[iParam]] %*% iIB.Lambda
        }
    }

    if("Psi" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$Psi)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            dOmega[[iParam]] <- tLambda.tiIB %*% skeleton$dmat.dparam$Psi[[iParam]] %*% iIB.Lambda
        }
    }

    if("sigma2" %in% names(skeleton$dmat.dparam)){
        dOmega[["sigma2"]] <- Sigma/param["sigma2"]
    }

    if("sigma2k" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$sigma2k)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            iFactor <- skeleton$dmat.dparam$sigma2k[[iParam]][,,"X"]+skeleton$dmat.dparam$sigma2k[[iParam]][,,"Y"]
            dOmega[[iParam]] <- Sigma*iFactor/param[iParam]
        }
    }

    if("cor" %in% names(skeleton$dmat.dparam)){
        iName.param <- names(skeleton$dmat.dparam$cor)
        for(iParam in iName.param){ ## iParam <- iName.param[1]
            if(abs(param[iParam])>1e-10){
                dOmega[[iParam]] <- Sigma*skeleton$dmat.dparam$cor[[iParam]]/param[iParam]
            }else{
                dOmega[[iParam]] <- Sigma*0
            }
        }
    }

    ## *** Export
    return(list(dmu = dmu, dOmega = dOmega))

}


## * updateD2Moment
#' @rdname updateD2Moment
updateD2Moment <- function(moment, skeleton, param){
    if(lava.options()$debug){cat("updateD2Moment \n")}

    ## ** Import information
    d2mu <- skeleton$d2mu.dparam
    d2Omega <- skeleton$d2Omega.dparam

    dalpha <- skeleton$dmat.dparam$alpha
    dLambda <- skeleton$dmat.dparam$Lambda
    dGamma <- skeleton$dmat.dparam$Gamma
    dB <- skeleton$dmat.dparam$B
    dPsi <- skeleton$dmat.dparam$Psi
    dsigma2 <- skeleton$dmat.dparam$sigma2
    dsigma2k <- skeleton$dmat.dparam$sigma2k
    dcor <- skeleton$dmat.dparam$cor

    Psi <- moment$Psi
    Lambda <- moment$Lambda
    iIB <- moment$iIB
    Psi.iIB <- moment$Psi.iIB
    iIB.Lambda <- moment$iIB.Lambda
    tLambda.tiIB <- moment$tLambda.tiIB
    alpha.XGamma.iIB <- moment$alpha.XGamma.iIB
    tiIB.Psi.iIB <- moment$tiIB.Psi.iIB
    tiIB.Psi.iIB.Lambda <- moment$tiIB.Psi.iIB.Lambda
    tLambda.tiIB.Psi.iIB <- moment$tLambda.tiIB.Psi.iIB
    Sigma <- moment$Sigma
    attr(Sigma,"detail") <- NULL
    SigmaValue <- attr(moment$Sigma,"detail")
    
    grid.mean <- skeleton$grid.d2moment$mean
    grid.var <- skeleton$grid.d2moment$var
    names.grid.mean <- names(grid.mean)
    names.grid.var <- names(grid.var)
    
    ## ** Compute partial derivative regarding the mean
    if("alpha.B" %in% names.grid.mean){
        for(iP in 1:NROW(grid.mean$alpha.B)){ # iP <- 1
            iName1 <- grid.mean$alpha.B[iP,"alpha"]
            iName2 <- grid.mean$alpha.B[iP,"B"]
            d2mu[[iName1]][[iName2]] <- dalpha[[iName1]] %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
        }
    }
    
    if("alpha.Lambda" %in% names.grid.mean){
        for(iP in 1:NROW(grid.mean$alpha.Lambda)){ # iP <- 1
            iName1 <- grid.mean$alpha.Lambda[iP,"alpha"]
            iName2 <- grid.mean$alpha.Lambda[iP,"Lambda"]
            d2mu[[iName1]][[iName2]] <- dalpha[[iName1]] %*% iIB %*% dLambda[[iName2]]
        }
    }

    if("Gamma.B" %in% names.grid.mean){
        for(iP in 1:NROW(grid.mean$Gamma.B)){ # iP <- 1
                iName1 <- grid.mean$Gamma.B[iP,"Gamma"]
                iName2 <- grid.mean$Gamma.B[iP,"B"]
                d2mu[[iName1]][[iName2]] <- dGamma[[iName1]] %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
        }
    }

    if("Gamma.Lambda" %in% names.grid.mean){
        for(iP in 1:NROW(grid.mean$Gamma.Lambda)){ # iP <- 1
                iName1 <- grid.mean$Gamma.Lambda[iP,"Gamma"]
                iName2 <- grid.mean$Gamma.Lambda[iP,"Lambda"]                
                d2mu[[iName1]][[iName2]] <- dGamma[[iName1]] %*% iIB %*% dLambda[[iName2]]
        }
    }

    if("Lambda.B" %in% names.grid.mean){
        for(iP in 1:NROW(grid.mean$Lambda.B)){ # iP <- 1
                iName1 <- grid.mean$Lambda.B[iP,"Lambda"]
                iName2 <- grid.mean$Lambda.B[iP,"B"]
                d2mu[[iName1]][[iName2]] <- alpha.XGamma.iIB %*% dB[[iName2]] %*% iIB %*% dLambda[[iName1]]
        }
    }

    if("B.B" %in% names.grid.mean){
        for(iP in 1:NROW(grid.mean$B.B)){ # iP <- 1
                iName1 <- grid.mean$B.B[iP,"B1"]
                iName2 <- grid.mean$B.B[iP,"B2"]

                term1 <- alpha.XGamma.iIB %*% dB[[iName2]] %*% iIB %*% dB[[iName1]] %*% iIB.Lambda
                term2 <- alpha.XGamma.iIB %*% dB[[iName1]] %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
                d2mu[[iName1]][[iName2]] <- term1 + term2
        }
    }

    ## ** Compute partial derivative regarding the variance
    if("Psi.Lambda" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$Psi.Lambda)){ # iP <- 1
            iName1 <- grid.var$Psi.Lambda[iP,"Psi"]
            iName2 <- grid.var$Psi.Lambda[iP,"Lambda"]

            term1 <- t(dLambda[[iName2]]) %*% t(iIB) %*% dPsi[[iName1]] %*% iIB.Lambda                
            d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
        }
    }
    
    if("Psi.B" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$Psi.B)){ # iP <- 1
                iName1 <- grid.var$Psi.B[iP,"Psi"]
                iName2 <- grid.var$Psi.B[iP,"B"]

                term1 <- t(iIB.Lambda) %*% t(dB[[iName2]]) %*% t(iIB) %*% dPsi[[iName1]] %*% iIB.Lambda
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
        }
    }
    
    if("Lambda.B" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$Lambda.B)){ # iP <- 1
                iName1 <- grid.var$Lambda.B[iP,"Lambda"]
                iName2 <- grid.var$Lambda.B[iP,"B"]

                term1 <- t(dLambda[[iName1]]) %*% t(iIB) %*% t(dB[[iName2]]) %*% t(iIB) %*% Psi %*% iIB.Lambda
                term2 <- t(dLambda[[iName1]]) %*% t(iIB) %*% Psi %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
                ## term2 <- tLambda.tiIB.Psi.iIB %*% dB[[iName2]] %*% iIB %*% dLambda[[iName1]]                
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1) + term2 + t(term2)
        }
    }

    if("Lambda.Lambda" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$Lambda.Lambda)){ # iP <- 1
            iName1 <- grid.var$Lambda.Lambda[iP,"Lambda1"]
            iName2 <- grid.var$Lambda.Lambda[iP,"Lambda2"]
                
            term1 <- t(dLambda[[iName1]]) %*% t(iIB) %*% Psi.iIB %*% dLambda[[iName2]]
            d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
        }
    }

    if("B.B" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$B.B)){ # iP <- 1
            iName1 <- grid.var$B.B[iP,"B1"]
            iName2 <- grid.var$B.B[iP,"B2"]

            term1 <- t(iIB.Lambda) %*% t(dB[[iName2]]) %*% t(iIB) %*% t(dB[[iName1]]) %*% t(iIB) %*% Psi.iIB %*% Lambda
            term2 <- t(iIB.Lambda) %*% t(dB[[iName1]]) %*% t(iIB) %*% t(dB[[iName2]]) %*% t(iIB) %*% Psi.iIB %*% Lambda
            term3 <- t(iIB.Lambda) %*% t(dB[[iName1]]) %*% t(iIB) %*% Psi.iIB %*% dB[[iName2]] %*% iIB %*% Lambda
            d2Omega[[iName1]][[iName2]] <- term1 + t(term1) + term2 + t(term2) + term3 + t(term3)
        }
    }

    if("sigma2.cor" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$sigma2.cor)){ # iP <- 1
            iName1 <- grid.var$sigma2.cor[iP,"sigma2"]
            iName2 <- grid.var$sigma2.cor[iP,"cor"]
            d2Omega[[iName1]][[iName2]] <- Sigma * (dsigma2[[iName1]]/param[iName1]) * (dcor[[iName2]]/param[iName2])
        }
    }

    if("sigma2.sigma2k" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$sigma2.sigma2k)){ # iP <- 1
            iName1 <- grid.var$sigma2.sigma2k[iP,"sigma2"]
            iName2 <- grid.var$sigma2.sigma2k[iP,"sigma2k"]

            iFactor <- dsigma2k[[iName2]][,,"X"] + dsigma2k[[iName2]][,,"Y"]
            d2Omega[[iName1]][[iName2]] <- Sigma * (dsigma2[[iName1]]/param[iName1]) * (iFactor/param[iName2])
        }
    }

    if("cor.sigma2k" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$cor.sigma2k)){ # iP <- 1
            iName1 <- grid.var$cor.sigma2k[iP,"cor"]
            iName2 <- grid.var$cor.sigma2k[iP,"sigma2k"]

            iFactor <- dsigma2k[[iName2]][,,"X"] + dsigma2k[[iName2]][,,"Y"]
            d2Omega[[iName1]][[iName2]] <- Sigma * (dcor[[iName1]]/param[iName1]) * (iFactor/param[iName2])
        }
    }

    if("sigma2k.sigma2k" %in% names.grid.var){
        for(iP in 1:NROW(grid.var$sigma2k.sigma2k)){ # iP <- 3
            iName1 <- grid.var$sigma2k.sigma2k[iP,"sigma2k1"]
            iName2 <- grid.var$sigma2k.sigma2k[iP,"sigma2k2"]
            iFactor <- dsigma2k[[iName1]][,,"X"] * dsigma2k[[iName2]][,,"Y"]
            if(iName1==iName2){
                iFactor <- 2 * iFactor
            }
            d2Omega[[iName1]][[iName2]] <- Sigma * (iFactor/(param[iName1]*param[iName2]))
        }
    }

    ## ** Export
    return(list(d2mu = d2mu, d2Omega = d2Omega))

}



######################################################################
### sCorrect-updateMoment.R ends here