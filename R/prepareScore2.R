### prepareScore2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: nov  9 2017 (16:23) 
##           By: Brice Ozenne
##     Update #: 445
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Prepare the computation of score2.
#' @description Compute partial derivatives regarding to the mean and the variance, and compute the design matrices.
#' @name prepareScore2
#' 
#' @param object a latent variable model
#' @param data [optional] data set.
#' @param name.endogenous [optional] name of the endogenous variables
#' @param name.latent [optional] name of the latent variables
#' 
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' e <- estimate(m,sim(m,1e2))
#' res <- prepareScore2(e)
#' res$skeleton$dt.param
#' @export
`prepareScore2` <-
  function(object, ...) UseMethod("prepareScore2")

## * prepareScore2.gls
#' @rdname prepareScore2
#' @export
prepareScore2.gls <- function(object, X, Omega,
                              var.coef, cor.coef,
                              n.cluster, n.endogenous, name.endogenous, index.obs){

### ** prepare
    class.var <- class(object$modelStruct$varStruct)
    class.cor <- class(object$modelStruct$corStruct)
    name.varcoef <- names(var.coef)
    name.corcoef <- names(cor.coef)
    n.varcoef <- length(var.coef)
    n.corcoef <- length(cor.coef)
    
    if("corSymm" %in% class.cor){
        M.corcoef <- matrix("", n.endogenous, n.endogenous,
                            dimnames = list(name.endogenous,name.endogenous))
        M.corcoef[lower.tri(M.corcoef)] <- name.corcoef
        M.corcoef[upper.tri(M.corcoef)] <- name.corcoef
    }
        
### ** score - mean
    name.X <- colnames(X)
    dmu.dtheta <- lapply(name.X, function(iCoef){
        M.tempo <- matrix(NA, nrow = n.cluster, ncol = n.endogenous)    
        M.tempo[index.obs] <- X[,iCoef]
        colnames(M.tempo) <- name.endogenous
        return(M.tempo)
    })
    names(dmu.dtheta) <- name.X

### ** score - variance/covariance
    dOmega.dtheta <- vector(mode = "list", length = n.corcoef + n.varcoef)
    names(dOmega.dtheta) <- c(name.corcoef, name.varcoef)

    for(iC in 1:n.cluster){ # iC <- 1
        iOmega <- Omega[[iC]]
        iSigma2.base <- diag(iOmega)
        iN.endogenous <- length(iSigma2.base)
        iName.endogenous <- colnames(iOmega)

        ## *** sigma2        
        dOmega.dtheta[["sigma2"]][[iC]] <- diag(iSigma2.base/var.coef["sigma2"],
                                                nrow = iN.endogenous, ncol = iN.endogenous)
        if("NULL" %in% class.cor == FALSE){
            dOmega.dtheta[["sigma2"]][[iC]][lower.tri(iOmega)] <- iOmega[lower.tri(iOmega)]/var.coef["sigma2"]
            dOmega.dtheta[["sigma2"]][[iC]][upper.tri(iOmega)] <- iOmega[upper.tri(iOmega)]/var.coef["sigma2"]
        }

        if(n.endogenous > 1){
            colnames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
            rownames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous 
        }
        
        ## *** other sigma
        if("varIdent" %in% class.var){
            for(iVar in setdiff(name.varcoef,"sigma2")){ # iVar <- name.varcoef[2]
                index.iVar <- iName.endogenous %in% iVar 
                dOmega.dtheta[[iVar]][[iC]] <- var.coef["sigma2"]*diag(index.iVar,iN.endogenous,iN.endogenous)

                if("NULL" %in% class.cor == FALSE){
                    index2.iVar <- setdiff(1:iN.endogenous,which(index.iVar))                
                    dOmega.dtheta[[iVar]][[iC]][index2.iVar,index.iVar] <- iOmega[index2.iVar,index.iVar]/var.coef[iVar]
                    dOmega.dtheta[[iVar]][[iC]][index.iVar,index2.iVar] <- iOmega[index.iVar,index2.iVar]/var.coef[iVar]
                }
                colnames(dOmega.dtheta[[iVar]][[iC]]) <- iName.endogenous
                rownames(dOmega.dtheta[[iVar]][[iC]]) <- iName.endogenous
            }
        }

        ## *** correlation coefficients
        if("corCompSymm" %in% class.cor){            
            dOmega.dtheta[[name.corcoef]][[iC]] <- sqrt(iSigma2.base) %*% t(sqrt(iSigma2.base)) - diag(diag(iOmega))
            colnames(dOmega.dtheta[[name.corcoef]][[iC]]) <- iName.endogenous
            rownames(dOmega.dtheta[[name.corcoef]][[iC]]) <- iName.endogenous
        }else if("corSymm" %in% class.cor){
            iM.corcoef <- M.corcoef[iName.endogenous,iName.endogenous,drop=FALSE]
            for(iVar in name.corcoef){ # iVar <- name.corcoef[1]
                dOmega.dtheta[[iVar]][[iC]] <- matrix(0, nrow = iN.endogenous, ncol = iN.endogenous)
                index.iVar <- which(iM.corcoef == iVar)
                dOmega.dtheta[[iVar]][[iC]][index.iVar] <- iOmega[index.iVar]/cor.coef[iVar]
                colnames(dOmega.dtheta[[iVar]][[iC]]) <- iName.endogenous
                rownames(dOmega.dtheta[[iVar]][[iC]]) <- iName.endogenous
            }
        }
    }
    
### ** export
    Omega_chol <- lapply(Omega,function(i){
        M <- chol(i)
        colnames(M) <- colnames(i)
        rownames(M) <- rownames(i)
        return(M)
    })
    iOmega <- lapply(Omega_chol,function(i){
        M <- chol2inv(i)
        colnames(M) <- colnames(i)
        rownames(M) <- rownames(i)
        return(M)
    }) ## faster compared to solve

    return(list(dmu.dtheta = dmu.dtheta,
                dOmega.dtheta = dOmega.dtheta,
                Omega_chol = Omega_chol,
                iOmega = iOmega))
    
}

## * prepareScore2.lme
#' @rdname prepareScore2
#' @export
prepareScore2.lme <- function(object, X, Omega,
                              var.coef, cor.coef,
                              n.cluster, n.endogenous, name.endogenous, index.obs){

### ** prepare
    class.var <- class(object$modelStruct$varStruct)
    name.varcoef <- names(var.coef)
    name.corcoef <- names(cor.coef)
    n.varcoef <- length(var.coef)
    n.corcoef <- length(cor.coef)
    
### ** score - mean
    name.X <- colnames(X)
    dmu.dtheta <- lapply(name.X, function(iCoef){
        M.tempo <- matrix(NA, nrow = n.cluster, ncol = n.endogenous)    
        M.tempo[index.obs] <- X[,iCoef]
        colnames(M.tempo) <- name.endogenous
        return(M.tempo)
    })
    names(dmu.dtheta) <- name.X

### ** score - variance/covariance
    dOmega.dtheta <- vector(mode = "list", length = n.corcoef + n.varcoef)
    names(dOmega.dtheta) <- c(name.corcoef, name.varcoef)

    for(iC in 1:n.cluster){ # iC <- 1
        iOmega <- Omega[[iC]]
        iSigma2.base <- diag(iOmega)-cor.coef
        iN.endogenous <- length(iSigma2.base)
        iName.endogenous <- colnames(iOmega)

        ## *** sigma2
        dOmega.dtheta[["sigma2"]][[iC]] <- diag(iSigma2.base/var.coef["sigma2"],iN.endogenous,iN.endogenous)        
        colnames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
        rownames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
        
        ## *** other sigma
        if("varIdent" %in% class.var){
            for(iVar in setdiff(name.varcoef,"sigma2")){ # iVar <- name.varcoef[2]
                index.iVar <- iName.endogenous %in% iVar 
                dOmega.dtheta[[iVar]][[iC]] <- var.coef["sigma2"]*diag(index.iVar,iN.endogenous,iN.endogenous)
                colnames(dOmega.dtheta[[iVar]][[iC]]) <- iName.endogenous
                rownames(dOmega.dtheta[[iVar]][[iC]]) <- iName.endogenous
            }
        }

        ## *** correlation coefficients
        dOmega.dtheta[[name.corcoef]][[iC]] <- matrix(1, nrow = iN.endogenous, ncol = iN.endogenous,
                                                      dimnames = list(iName.endogenous,iName.endogenous))
    }
    
### ** export
    Omega_chol <- lapply(Omega,function(i){
        M <- chol(i)
        colnames(M) <- colnames(i)
        rownames(M) <- rownames(i)
        return(M)
    })
    iOmega <- lapply(Omega_chol,function(i){
        M <- chol2inv(i)
        colnames(M) <- colnames(i)
        rownames(M) <- rownames(i)
        return(M)
    }) ## faster compared to solve

    return(list(dmu.dtheta = dmu.dtheta,
                dOmega.dtheta = dOmega.dtheta,
                Omega_chol = Omega_chol,
                iOmega = iOmega))
    
}
    



#----------------------------------------------------------------------
### prepareScore2.R ends here


#' @rdname prepareScore2
#' @export
`prepareScore2<-` <-
  function(object, ..., value) UseMethod("prepareScore2<-")

#' @rdname prepareScore2
#' @export
"prepareScore2<-.lvmfit" <- function(object, ..., value) {
    object$prepareScore2  <- prepareScore2(lava::Model(object), data = value, ...)
    return(object)
}
## * prepareScore2.lvm
#' @rdname prepareScore2
#' @export
prepareScore2.lvm <- function(object, data,
                              name.endogenous = NULL, name.latent = NULL){

### ** normalize arguments
    if(is.null(name.endogenous)){name.endogenous <- endogenous(object)}
    n.endogenous <- length(name.endogenous)
    if(is.null(name.latent)){name.latent <- latent(object)}
    n.latent <- length(name.latent)

    if(!is.matrix(data)){
        data <- as.matrix(data)
    }

    prepareScore2 <- list()
    
### ** Compute skeleton   
    prepareScore2$skeleton <- skeleton(object,
                                       name.endogenous = name.endogenous, 
                                       name.latent = name.latent, 
                                       as.lava = TRUE)
    
### ** Initialize partial derivatives
    prepareScore2$dtheta <- skeletonDtheta(object, data = data,
                                           dt.param.all = prepareScore2$skeleton$dt.param,
                                           param2originalLink = prepareScore2$skeleton$param2originalLink,
                                           name.endogenous = name.endogenous, 
                                           name.latent = name.latent)
    
### ** Export
    return(prepareScore2)
}
    
    
## * prepareScore2.lvmfit
#' @rdname prepareScore2
#' @export
prepareScore2.lvmfit <- function(object, data = NULL, p = NULL,
                                 name.endogenous = NULL, name.latent = NULL){

### ** normalize arguments
    if(is.null(name.endogenous)){name.endogenous <- endogenous(object)}
    n.endogenous <- length(name.endogenous)
    if(is.null(name.latent)){name.latent <- latent(object)}
    n.latent <- length(name.latent)

    if(is.null(data)){
        data <- model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }

    if(is.null(p)){
        p <- pars(object)        
    }
    
    prepareScore2 <- list()

### ** Update skeleton with current estimates
    prepareScore2$skeleton <- skeleton(object, data = data, p = p,
                                       name.endogenous = name.endogenous, 
                                       name.latent = name.latent, 
                                       as.lava = TRUE)
    
### ** Update partial derivatives with current estimates
    prepareScore2$dtheta <- skeletonDtheta(object, data = data,
                                           dt.param.all = prepareScore2$skeleton$dt.param,
                                           param2originalLink = prepareScore2$skeleton$param2originalLink,
                                           name.endogenous = name.endogenous, 
                                           name.latent = name.latent,
                                           B = prepareScore2$skeleton$value$B,
                                           alpha.XGamma = prepareScore2$skeleton$value$alpha.XGamma,
                                           Lambda = prepareScore2$skeleton$value$Lambda,
                                           Psi = prepareScore2$skeleton$value$Psi)
    
### ** Export
    return(prepareScore2)    
    
}

