### sCorrect-updateMoment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 10 2019 (09:58) 
## Version: 
## Last-Updated: dec 10 2019 (17:29) 
##           By: Brice Ozenne
##     Update #: 44
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
                         param, endogenous, latent, n.cluster){
    if(TRUE){cat("updateMoment \n")}


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
            }else{
                index.update <- which(!is.na(skeleton[[iUpdate]]))
                value[[iUpdate]][index.update] <- param[skeleton[[iUpdate]][index.update]]
            }
        }
    }

    ## ** Pre-compute relevant quantities
    if(!is.null(value$B)){
        value$iIB <- solve(diag(1,n.latent,n.latent)-value$B)
    }else{
        value$iIB <- diag(1,n.latent,n.latent)
    }
    if(!is.null(value$Lambda)){
        value$iIB.Lambda <-  value$iIB %*% value$Lambda    
    }
    if(!is.null(value$Psi)){
        value$Psi.iIB <- value$Psi %*% value$iIB
        value$tLambda.tiIB.Psi.iIB <- t(value$iIB.Lambda) %*% value$Psi.iIB
    }

    ## ** Compute mean
    value$mu <- matrix(0, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL,endogenous))
    if(!is.null(value$nu)){
        value$mu <- value$mu + sweep(skeleton$Xnu, MARGIN = 2, FUN = "*", STATS = value$nu)
    }
    if(!is.null(value$K)){
        value$mu <- value$mu + do.call(cbind,lapply(endogenous, function(iE){skeleton$XK[[iE]] %*% value$K[,iE]})) ## iE <- endogenous[1]
    }
    if(!is.null(value$alpha)){
        value$mu <- value$mu + skeleton$Xalpha %o% (value$alpha %*% value$iIB.Lambda)[1,]
    }
    if(!is.null(value$Gamma)){
        value$mu <- value$mu + do.call(cbind,lapply(latent, function(iL){skeleton$XGamma[[iL]] %*% value$Gamma[,iL]})) %*% value$iIB.Lambda
    }

    ## ** Compute variance
    value$Omega <- matrix(0, nrow = n.endogenous, ncol = n.endogenous, 
                          dimnames = list(endogenous,endogenous))
    
    if(!is.null(value$Sigma)){
        value$Omega <- value$Omega + value$Sigma
    }
    if(!is.null(value$Psi)){
        value$Omega <- value$Omega + value$tLambda.tiIB.Psi.iIB %*% value$Lambda
    }
    ## ** Export
    return(value)
}

## * updateDMoment
#' @rdname updateMoment
updateDMoment <- function(moment, skeleton,
                          toUpdate,
                          endogenous, latent, n.cluster){

    n.endogenous <- length(endogenous)
    n.latent <- length(latent)

    ## ** Compute partial derivative regarding for each matrix of coefficients
    skeleton$dmu.dparam

    skeleton$dmat.dparam
    
    moment
    browser()
    
    ## ** Compute partial derivative regarding the mean
    
    ## ** Compute partial derivative regarding the variance
    
    ## from Moment
    type <- object$conditionalMoment$skeleton$type
    iIB.Lambda <- object$conditionalMoment$value$iIB.Lambda
    alpha.XGamma.iIB <- object$conditionalMoment$value$alpha.XGamma.iIB
    tLambda.tiIB.Psi.iIB <- object$conditionalMoment$value$tLambda.tiIB.Psi.iIB

    ## from dMoment.init
    dmu <- object$conditionalMoment$dMoment.init$dmu
    dOmega <- object$conditionalMoment$dMoment.init$dOmega
    dLambda <- object$conditionalMoment$dMoment.init$dLambda
    dB <- object$conditionalMoment$dMoment.init$dB
    dPsi <- object$conditionalMoment$dMoment.init$dPsi
    toUpdate <- object$conditionalMoment$dMoment.init$toUpdate
    name2Update <- names(toUpdate)
    type2Update <- type[name2Update]
    
    ## *** Update partial derivatives

    ## **** mean coefficients
    type2Update.meanparam <- type2Update[type2Update %in% c("alpha","Lambda","Gamma","B")]
    name2Update.meanparam <- names(type2Update.meanparam)
    n2Update.meanparam <- length(name2Update.meanparam)
        
    if(n2Update.meanparam>0){
        for(iP in 1:n2Update.meanparam){ # iP <- 1
            iType <- type2Update.meanparam[iP]
            iName <- name2Update.meanparam[iP]
            
            if(iType == "alpha"){
                dmu[[iName]] <- dmu[[iName]] %*% iIB.Lambda
            }else if(iType == "Gamma"){
                dmu[[iName]] <- dmu[[iName]] %*% iIB.Lambda 
            }else if(iType == "Lambda"){
                dmu[[iName]] <- alpha.XGamma.iIB %*% dLambda[[iName]]
            }else if(iType == "B"){
                dmu[[iName]] <- alpha.XGamma.iIB %*% dB[[iName]] %*% iIB.Lambda
            }

            colnames(dmu[[iName]]) <- name.endogenous
        }
    }

    ## **** variance-covariance coefficients
    type2Update.vcovparam <- type2Update[type2Update %in% c("Psi_var","Psi_cov","Lambda","B")]
    name2Update.vcovparam <- names(type2Update.vcovparam)
    n2Update.vcovparam <- length(name2Update.vcovparam)

    if(n2Update.vcovparam>0){
        for(iP in 1:n2Update.vcovparam){ # iP <- 1
            iType <- type2Update.vcovparam[iP]
            iName <- name2Update.vcovparam[iP]
        
            if(iType %in% "Psi_var"){
                dOmega[[iName]] <-  t(iIB.Lambda) %*% dPsi[[iName]] %*% iIB.Lambda
            }else if(iType %in% "Psi_cov"){
                dOmega[[iName]] <-  t(iIB.Lambda) %*% dPsi[[iName]] %*% iIB.Lambda
            }else if(iType == "Lambda"){
                dOmega[[iName]] <- tLambda.tiIB.Psi.iIB %*% dLambda[[iName]]
                dOmega[[iName]] <- dOmega[[iName]] + t(dOmega[[iName]])
            }else if(iType == "B"){
                dOmega[[iName]] <- tLambda.tiIB.Psi.iIB %*% dB[[iName]] %*% iIB.Lambda
                dOmega[[iName]] <- dOmega[[iName]] + t(dOmega[[iName]])
            }

            colnames(dOmega[[iName]]) <- name.endogenous
            rownames(dOmega[[iName]]) <- name.endogenous
        }
    }

    ## *** Export
    return(list(dmu = dmu, dOmega = dOmega))

}


## * updateD2Moment
#' @rdname updateD2Moment
updateD2Moment <- function(){
    
    ## *** Import information
    n.endogenous <- NCOL(object$conditionalMoment$Omega)

    ## from Moment
    Psi <- object$conditionalMoment$value$Psi
    Lambda <- object$conditionalMoment$value$Lambda
    iIB <- object$conditionalMoment$value$iIB
    Psi.iIB <- object$conditionalMoment$value$Psi.iIB
    iIB.Lambda <- object$conditionalMoment$value$iIB.Lambda
    alpha.XGamma.iIB <- object$conditionalMoment$value$alpha.XGamma.iIB
    type <- object$conditionalMoment$skeleton$type

    ## from dMoment.init
    dLambda <- object$conditionalMoment$dMoment.init$dLambda
    dB <- object$conditionalMoment$dMoment.init$dB
    dPsi <- object$conditionalMoment$dMoment.init$dPsi
    
    ## from d2Moment.init
    d2mu <- object$conditionalMoment$d2Moment.init$d2mu
    d2Omega <- object$conditionalMoment$d2Moment.init$d2Omega

    grid.mean <- object$conditionalMoment$d2Moment.init$grid.mean
    grid.vcov <- object$conditionalMoment$d2Moment.init$grid.vcov

    n.mean <- object$conditionalMoment$d2Moment.init$n.mean
    n.vcov <- object$conditionalMoment$d2Moment.init$n.vcov

    toUpdate <- object$conditionalMoment$d2Moment.init$toUpdate
    ##    names(object$conditionalMoment$d2Moment)
    
    ## *** second order partial derivatives
    if(any(toUpdate)){
        
        ## **** mean coefficients        
        if(toUpdate["alpha.B"]){
            for(iP in 1:n.mean$alpha.B){ # iP <- 1
                iName1 <- grid.mean$alpha.B[iP,"alpha"]
                iName2 <- grid.mean$alpha.B[iP,"B"]

                d2mu[[iName1]][[iName2]] <- d2mu[[iName1]][[iName2]] %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
            }
        }
        
        if(toUpdate["alpha.Lambda"]){
            for(iP in 1:n.mean$alpha.Lambda){ # iP <- 1
                iName1 <- grid.mean$alpha.Lambda[iP,"alpha"]
                iName2 <- grid.mean$alpha.Lambda[iP,"Lambda"]

                d2mu[[iName1]][[iName2]] <- d2mu[[iName1]][[iName2]] %*% iIB %*% dLambda[[iName2]]
                
            }
        }

        if(toUpdate["Gamma.B"]){
            for(iP in 1:n.mean$Gamma.B){ # iP <- 1
                iName1 <- grid.mean$Gamma.B[iP,"Gamma"]
                iName2 <- grid.mean$Gamma.B[iP,"B"]

                d2mu[[iName1]][[iName2]] <- d2mu[[iName1]][[iName2]] %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
            }
        }        

        if(toUpdate["Gamma.Lambda"]){
            for(iP in 1:n.mean$Gamma.Lambda){ # iP <- 1
                iName1 <- grid.mean$Gamma.Lambda[iP,"Gamma"]
                iName2 <- grid.mean$Gamma.Lambda[iP,"Lambda"]                

                d2mu[[iName1]][[iName2]] <- d2mu[[iName1]][[iName2]] %*% iIB %*% dLambda[[iName2]]
            }
        }        

        if(toUpdate["Lambda.B"]){
            for(iP in 1:n.mean$Lambda.B){ # iP <- 1
                iName1 <- grid.mean$Lambda.B[iP,"Lambda"]
                iName2 <- grid.mean$Lambda.B[iP,"B"]

                d2mu[[iName1]][[iName2]] <- alpha.XGamma.iIB %*% dB[[iName2]] %*% iIB %*% dLambda[[iName1]]
            }
        }

        if(toUpdate["B.B"]){
            for(iP in 1:n.mean$B.B){ # iP <- 1
                iName1 <- grid.mean$B.B[iP,"B1"]
                iName2 <- grid.mean$B.B[iP,"B2"]

                term1 <- alpha.XGamma.iIB %*% dB[[iName2]] %*% iIB %*% dB[[iName1]] %*% iIB.Lambda
                term2 <- alpha.XGamma.iIB %*% dB[[iName1]] %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
                d2mu[[iName1]][[iName2]] <- term1 + term2
            }
        }

        ## **** variance-covariance coefficients
        if(toUpdate["Psi.Lambda"]){
            for(iP in 1:n.vcov$Psi.Lambda){ # iP <- 1
                iName1 <- grid.vcov$Psi.Lambda[iP,"Psi"]
                iName2 <- grid.vcov$Psi.Lambda[iP,"Lambda"]

                term1 <- t(dLambda[[iName2]]) %*% t(iIB) %*% dPsi[[iName1]] %*% iIB.Lambda                
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
            }
        }

        if(toUpdate["Psi.B"]){
            for(iP in 1:n.vcov$Psi.B){ # iP <- 1
                iName1 <- grid.vcov$Psi.B[iP,"Psi"]
                iName2 <- grid.vcov$Psi.B[iP,"B"]

                term1 <- t(iIB.Lambda) %*% t(dB[[iName2]]) %*% t(iIB) %*% dPsi[[iName1]] %*% iIB.Lambda
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
            }
        }

        if(toUpdate["Lambda.B"]){
            for(iP in 1:n.vcov$Lambda.B){ # iP <- 1
                iName1 <- grid.vcov$Lambda.B[iP,"Lambda"]
                iName2 <- grid.vcov$Lambda.B[iP,"B"]

                term1 <- t(dLambda[[iName1]]) %*% t(iIB) %*% t(dB[[iName2]]) %*% t(iIB) %*% Psi %*% iIB.Lambda
                term2 <- t(dLambda[[iName1]]) %*% t(iIB) %*% Psi %*% iIB %*% dB[[iName2]] %*% iIB.Lambda
                ## term2 <- tLambda.tiIB.Psi.iIB %*% dB[[iName2]] %*% iIB %*% dLambda[[iName1]]                
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1) + term2 + t(term2)
            }
        }

        if(toUpdate["Lambda.Lambda"]){
            for(iP in 1:n.vcov$Lambda.Lambda){ # iP <- 1
                iName1 <- grid.vcov$Lambda.Lambda[iP,"Lambda1"]
                iName2 <- grid.vcov$Lambda.Lambda[iP,"Lambda2"]
                
                term1 <- t(dLambda[[iName1]]) %*% t(iIB) %*% Psi.iIB %*% dLambda[[iName2]]
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1)
            }
        }

        if(toUpdate["B.B"]){
            for(iP in 1:n.vcov$B.B){ # iP <- 1
                iName1 <- grid.vcov$B.B[iP,"B1"]
                iName2 <- grid.vcov$B.B[iP,"B2"]

                term1 <- t(iIB.Lambda) %*% t(dB[[iName2]]) %*% t(iIB) %*% t(dB[[iName1]]) %*% t(iIB) %*% Psi.iIB %*% Lambda
                term2 <- t(iIB.Lambda) %*% t(dB[[iName1]]) %*% t(iIB) %*% t(dB[[iName2]]) %*% t(iIB) %*% Psi.iIB %*% Lambda
                term3 <- t(iIB.Lambda) %*% t(dB[[iName1]]) %*% t(iIB) %*% Psi.iIB %*% dB[[iName2]] %*% iIB %*% Lambda
                d2Omega[[iName1]][[iName2]] <- term1 + t(term1) + term2 + t(term2) + term3 + t(term3)
            }
        }

    }

    ## *** Export
    return(list(d2mu = d2mu, d2Omega = d2Omega))

}



######################################################################
### sCorrect-updateMoment.R ends here
