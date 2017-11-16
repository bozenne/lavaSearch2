### prepareScore2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: nov 16 2017 (15:55) 
##           By: Brice Ozenne
##     Update #: 592
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
                              param, attr.param,
                              n.cluster, n.endogenous, name.endogenous, index.obs){
    
### ** prepare
    name.varcoef <- attr.param$var.coef
    name.corcoef <- attr.param$cor.coef
    n.varcoef <- length(name.varcoef)
    n.corcoef <- length(name.corcoef)
    var.coef <- param[name.varcoef]
    cor.coef <- param[name.corcoef]

    class.var <- class(object$modelStruct$varStruct)
    class.cor <- class(object$modelStruct$corStruct)
    
    ## *** diagonal terms
    name.other <- setdiff(names(var.coef),"sigma2")
    if("NULL" %in% class.var){            
        sigma2.base <- setNames(rep(var.coef["sigma2"],n.endogenous), name.endogenous)
    }else{
        sigma2.base <- setNames(var.coef["sigma2"]*c(1,var.coef[name.other]), name.endogenous)            
    }
    sigma2.base0 <- sigma2.base/var.coef["sigma2"]

    ## *** variance
    dOmega.dsigma2 <- diag(sigma2.base0[name.endogenous], nrow = n.endogenous, ncol = n.endogenous)
    dimnames(dOmega.dsigma2) <-  list(name.endogenous, name.endogenous)

    if("NULL" %in% class.var == FALSE){           
        dOmega.dsigmaOther <- lapply(name.other, function(x){
            M <- var.coef["sigma2"]*diag(name.endogenous %in% x, nrow = n.endogenous, ncol = n.endogenous)
            dimnames(M) <- list(name.endogenous, name.endogenous)
            return(M)
        })
        names(dOmega.dsigmaOther) <- name.other
    }
    
    ## *** correlation
    if("NULL" %in% class.cor == FALSE){
        M.corcoef <- matrix("", n.endogenous, n.endogenous,
                            dimnames = list(name.endogenous,name.endogenous))
        M.corcoef[which(lower.tri(M.corcoef))] <- name.corcoef
        M.corcoef <- symmetrize(M.corcoef)

        index.lower.tri <- which(lower.tri(M.corcoef))
        indexArr.lower.tri <- which(lower.tri(M.corcoef), arr.ind = TRUE)

        Msigma2.base0 <- matrix(0, n.endogenous, n.endogenous,
                                dimnames = list(name.endogenous, name.endogenous))
        Msigma2.base0[index.lower.tri] <- apply(indexArr.lower.tri, 1, function(x){sqrt(prod(sigma2.base0[x]))})
        Msigma2.base0 <- symmetrize(Msigma2.base0)
        
        ## dOmega.dsigma2
        dOmega.dsigma2[index.lower.tri] <- Msigma2.base0[index.lower.tri] * cor.coef[M.corcoef[index.lower.tri]]
        dOmega.dsigma2 <- symmetrize(dOmega.dsigma2)

        ## dOmega.dcor
        dOmega.dcor <- vector("list", length = n.corcoef)
        names(dOmega.dcor) <- name.corcoef
        for(iVar in name.corcoef){
            dOmega.dcor[[iVar]] <- Msigma2.base0 * var.coef["sigma2"] * (M.corcoef==iVar)
        }

        ## dOmega.dsigmaOther
        if("NULL" %in% class.var == FALSE){
            for(iVar in name.other){ # iVar <- name.other[1]
                iEndogenous <- which(name.endogenous==iVar)
                index.iVar <- which(rowSums(indexArr.lower.tri==iEndogenous)>0)

                ##  d sqrt(x) / d x = 1/(2 sqrt(x)) = sqrt(x) / (2*x)
                dOmega.dsigmaOther[[iVar]][index.lower.tri[index.iVar]] <- var.coef["sigma2"]*dOmega.dsigma2[index.lower.tri[index.iVar]]/(2*var.coef[iVar])
                dOmega.dsigmaOther[[iVar]] <- symmetrize(dOmega.dsigmaOther[[iVar]])
            
            }
        }
        
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
       
        iName.endogenous <- colnames(Omega[[iC]])
        iIndex.lower.tri <- lower.tri(Omega[[iC]])
        iN.endogenous <- length(iName.endogenous)

        ## *** sigma2
        dOmega.dtheta[["sigma2"]][[iC]] <- dOmega.dsigma2[iName.endogenous,iName.endogenous]
        
        ## *** other sigma
        if("NULL" %in% class.var == FALSE){
            for(iVar in name.other){ # iVar <- name.varcoef[2]
                dOmega.dtheta[[iVar]][[iC]] <- dOmega.dsigmaOther[[iVar]][iName.endogenous,iName.endogenous]
            }
        }

        ## *** correlation coefficients
        if("NULL" %in% class.cor == FALSE){           
            for(iVar in name.corcoef){ # iVar <- name.corcoef[1]
                dOmega.dtheta[[iVar]][[iC]] <- dOmega.dcor[[iVar]][iName.endogenous,iName.endogenous]
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
                              param, attr.param,
                              n.cluster, n.endogenous, name.endogenous, index.obs){

    resGLS <- prepareScore2.gls(object, X = X, Omega = Omega, param = param, attr.param = attr.param,
                                n.cluster = n.cluster,
                                n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                index.obs = index.obs)
        
### ** random effect
    name.rancoef <- attr.param$ran.coef
    resGLS$dOmega.dtheta[[name.rancoef]] <- vector("list",length=n.cluster)
    for(iC in 1:n.cluster){ # iC <- 1
        iName.endogenous <- colnames(Omega[[iC]])
        iN.endogenous <- length(iName.endogenous)
        resGLS$dOmega.dtheta[[name.rancoef]][[iC]] <- matrix(1, nrow = iN.endogenous, ncol = iN.endogenous,
                                                             dimnames = list(iName.endogenous,iName.endogenous))
    }


### ** export
    return(resGLS)
}
    

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




#----------------------------------------------------------------------
### prepareScore2.R ends here
