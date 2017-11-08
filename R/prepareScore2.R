### prepareScore2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: nov  8 2017 (21:56) 
##           By: Brice Ozenne
##     Update #: 356
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Prepare the computation of score2.
#' @description Prepare the computation of score2.
#' @name prepareScore2
#' 
#' @param object a latent variable model
#' @param data [optional] data set.
#' @param name.endogenous [optional] name of the endogenous variables
#' @param n.endogenous [optional] number of endogenous variables
#' @param name.latent [optional] name of the latent variables
#' @param n.latent [optional] number of latent variables
#' 
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' e <- estimate(m,sim(m,1e2))
#' prepareScore2(e)
#' @export
`prepareScore2` <-
  function(object, ...) UseMethod("prepareScore2")

## * prepareScore2.lvmfit
#' @rdname prepareScore2
#' @export
prepareScore2.lvmfit <- function(object, data = NULL,
                                 name.endogenous = NULL, n.endogenous = NULL,
                                 name.latent = NULL, n.latent = NULL){

### ** normalize arguments
    if(is.null(name.endogenous)){name.endogenous <- endogenous(object)}
    if(is.null(n.endogenous)){n.endogenous <- length(name.endogenous)}
    if(is.null(name.latent)){name.latent <- latent(object)}
    if(is.null(n.latent)){n.latent <- length(name.latent)}

    if(is.null(data)){
        data <- model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    
### ** Compute skeleton
    resSkeleton <- skeleton(object,
                            name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                            name.latent = name.latent, n.latent = n.latent,
                            as.lava = TRUE, update.value = FALSE)

    dt.param.all <- resSkeleton$dt.param
    dt.param <- dt.param.all[is.na(value) & factice == FALSE]
    Utype.by.detail <- dt.param[,.(n.type = length(unique(detail))), by = param][["n.type"]]
    if(any(Utype.by.detail>1)){
        stop("cannot constrain two parameters of different types to be equal \n")
    }
    name.param <- dt.param[!duplicated(param),param]
    n.param <- length(name.param)

    param2originalLink <- resSkeleton$param2originalLink
    name.originalLink <- as.character(param2originalLink)
    
### ** prepare
    n <- NROW(data)
    name.data <- colnames(data)
    
    mean.param <- c("nu","K","alpha","Gamma","Lambda","B")
    vcov.param <- c("Sigma_var","Sigma_cov","Psi_var","Psi_cov","Lambda","B")    
    dmu.dtheta <- list()
    dOmega.dtheta <- list()
    dLambda.dtheta <- list()
    dB.dtheta <- list()
    dPsi.dtheta <- list()

    type <- setNames(vector(mode = "character", n.param),name.originalLink)
    toUpdate <- setNames(vector(mode = "logical", n.param),name.originalLink)
    
### ** Compute derivative or prepare for the derivative
    for(iName in name.param){ # iName <- name.param[1]
        iName2 <- as.character(param2originalLink[iName])
        type[iName2] <- unique(dt.param[param == iName,detail])
        iY <- dt.param[param %in% iName,Y]
        iX <- dt.param[param %in% iName,X]

        ## *** derivative regarding the mean        
        if(type[iName2] %in% mean.param){
            if(type[iName2]=="nu"){
                dmu.dtheta[[iName2]] <- matrix(as.numeric(name.endogenous %in% iY),
                                              nrow = n, ncol = n.endogenous, byrow = TRUE,
                                              dimnames = list(NULL, name.endogenous))
                toUpdate[iName2] <- FALSE
            }else if(type[iName2]=="K"){
                dmu.dtheta[[iName2]] <- matrix(0, nrow = n, ncol = n.endogenous, byrow = TRUE,
                                              dimnames = list(NULL, name.endogenous))
                for(Y.tempo in unique(iY)){
                    dmu.dtheta[[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
                }
                toUpdate[iName2] <- FALSE
            }else if(type[iName2]=="alpha"){
                dmu.dtheta[[iName2]] <- matrix(as.numeric(name.latent %in% unique(iY)), nrow = n, ncol = n.latent, byrow = TRUE,
                                              dimnames = list(NULL, name.latent))                
                toUpdate[iName2] <- TRUE
            }else if(type[iName2]=="Gamma"){
                dmu.dtheta[[iName2]] <- matrix(0, nrow = n, ncol = n.latent, byrow = TRUE,
                                              dimnames = list(NULL, name.latent))
                for(Y.tempo in unique(iY)){ # Y.tempo <- "eta"
                    dmu.dtheta[[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
                }
                toUpdate[iName2] <- TRUE
            }
        }
        
        ## *** derivative regarding the residual variance covariance
        if(type[iName2] %in% vcov.param){
            
            if(type[iName2]=="Sigma_var"){
                dOmega.dtheta[[iName2]] <- matrix(0,
                                                  nrow = n.endogenous, ncol = n.endogenous, byrow = TRUE,
                                                  dimnames = list(name.endogenous, name.endogenous))
                dOmega.dtheta[[iName2]][match(iX, name.endogenous) + (match(iY, name.endogenous) - 1) * n.endogenous] <- 1
                toUpdate[iName2] <- FALSE
            }else if(type[iName2]=="Sigma_cov"){
                dOmega.dtheta[[iName2]] <- matrix(0,
                                                  nrow = n.endogenous, ncol = n.endogenous, byrow = TRUE,
                                                  dimnames = list(name.endogenous, name.endogenous))
                dOmega.dtheta[[iName2]][match(iX, name.endogenous) + (match(iY, name.endogenous) - 1) * n.endogenous] <- 1
                dOmega.dtheta[[iName2]][match(iY, name.endogenous) + (match(iX, name.endogenous) - 1) * n.endogenous] <- 1
                toUpdate[iName2] <- FALSE
            }
            
        }        

        ## *** matrices
        if(type[iName2]=="Lambda"){            
            dLambda.dtheta[[iName2]] <- matrix(0,
                                               nrow = n.latent, ncol = n.endogenous, byrow = TRUE,
                                               dimnames = list(name.latent, name.endogenous))
            dLambda.dtheta[[iName2]][match(iX, name.latent) + (match(iY, name.endogenous) - 1) * n.latent] <- 1            
            toUpdate[iName2] <- TRUE
        }else if(type[iName2]=="B"){
            dB.dtheta[[iName2]] <- matrix(0,
                                          nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                          dimnames = list(name.latent, name.latent))
            dB.dtheta[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            toUpdate[iName2] <- TRUE
        }else if(type[iName2]=="Psi_var"){
            dPsi.dtheta[[iName2]] <- matrix(0,
                                            nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                            dimnames = list(name.latent, name.latent))
            dPsi.dtheta[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            toUpdate[iName2] <- TRUE
        }else if(type[iName2]=="Psi_cov"){
            dPsi.dtheta[[iName2]] <- matrix(0,
                                            nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                            dimnames = list(name.latent, name.latent))
            dPsi.dtheta[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            dPsi.dtheta[[iName2]][match(iY, name.latent) + (match(iX, name.latent) - 1) * n.latent] <- 1            
            toUpdate[iName2] <- TRUE
        } 
    }
    
### ** export

    return(list(
        dmu.dtheta = dmu.dtheta,
        dOmega.dtheta = dOmega.dtheta,
        dLambda.dtheta = dLambda.dtheta,
        dB.dtheta = dB.dtheta,
        dPsi.dtheta = dPsi.dtheta,
        skeleton = resSkeleton$skeleton,
        value = resSkeleton$value,
        type = type,
        toUpdate = toUpdate
    ))
    
}

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
        dOmega.dtheta[["sigma2"]][[iC]] <- diag(iSigma2.base/var.coef["sigma2"],iN.endogenous,iN.endogenous)
        if("NULL" %in% class.cor){
            dOmega.dtheta[["sigma2"]][[iC]][lower.tri(iOmega)] <- iOmega[lower.tri(iOmega)]/var.coef["sigma2"]
            dOmega.dtheta[["sigma2"]][[iC]][upper.tri(iOmega)] <- iOmega[upper.tri(iOmega)]/var.coef["sigma2"]
            colnames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
            rownames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
        }
  
        ## *** other sigma
        if("NULL" %in% class.var){
            for(iVar in setdiff(name.varcoef,"sigma2")){ # iVar <- name.varcoef[2]
                index.iVar <- iName.endogenous %in% iVar 
                dOmega.dtheta[[iVar]][[iC]] <- var.coef["sigma2"]*diag(index.iVar,iN.endogenous,iN.endogenous)

                if("NULL" %in% class.cor){
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
        iSigma2.base <- diag(iOmega)
        iN.endogenous <- length(iSigma2.base)
        iName.endogenous <- colnames(iOmega)

        ## *** sigma2
        dOmega.dtheta[["sigma2"]][[iC]] <- diag(iSigma2.base/var.coef["sigma2"],iN.endogenous,iN.endogenous)
        dOmega.dtheta[["sigma2"]][[iC]][lower.tri(iOmega)] <- iOmega[lower.tri(iOmega)]/var.coef["sigma2"]
        dOmega.dtheta[["sigma2"]][[iC]][upper.tri(iOmega)] <- iOmega[upper.tri(iOmega)]/var.coef["sigma2"]
        colnames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
        rownames(dOmega.dtheta[["sigma2"]][[iC]]) <- iName.endogenous
        
        ## *** other sigma
        if("NULL" %in% class.var){
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
