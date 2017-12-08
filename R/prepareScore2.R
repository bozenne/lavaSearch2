### prepareScore2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: nov 30 2017 (16:09) 
##           By: Brice Ozenne
##     Update #: 639
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
prepareScore2.gls <- function(object, X, 
                              param, attr.param,
                              second.order,
                              n.cluster, n.endogenous, name.endogenous, index.obs){
    
### ** prepare
    
    ## *** coefficients
    name.varcoef <- attr.param$var.coef
    name.corcoef <- attr.param$cor.coef
    n.varcoef <- length(name.varcoef)
    n.corcoef <- length(name.corcoef)
    var.coef <- param[name.varcoef]
    cor.coef <- param[name.corcoef]

    class.var <- class(object$modelStruct$varStruct)
    class.cor <- class(object$modelStruct$corStruct)

    ## *** variance terms
    name.other <- setdiff(names(var.coef),"sigma2")
    if("NULL" %in% class.var){            
        sigma2.base <- setNames(rep(var.coef["sigma2"],n.endogenous), name.endogenous)
    }else{
        sigma2.base <- setNames(var.coef["sigma2"]*c(1,var.coef[name.other]), name.endogenous)            
    }
    sigma2.base0 <- sigma2.base/var.coef["sigma2"]

    ## *** corelation terms
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
    
    ## *** dispersion coefficient
    dOmega.dtheta[["sigma2"]] <- diag(sigma2.base0[name.endogenous], nrow = n.endogenous, ncol = n.endogenous)
   
    if("NULL" %in% class.cor == FALSE){
        dOmega.dtheta[["sigma2"]][index.lower.tri] <- Msigma2.base0[index.lower.tri] * cor.coef[M.corcoef[index.lower.tri]]
        dOmega.dtheta[["sigma2"]] <- symmetrize(dOmega.dtheta[["sigma2"]])      
    }
    dimnames(dOmega.dtheta[["sigma2"]]) <-  list(name.endogenous, name.endogenous)

    ## *** other variance coefficients
    if("NULL" %in% class.var == FALSE){
        for(iVar in name.other){
            dOmega.dtheta[[iVar]] <- var.coef["sigma2"]*diag(name.endogenous %in% iVar,
                                                             nrow = n.endogenous, ncol = n.endogenous)

            if("NULL" %in% class.cor == FALSE){
                iEndogenous <- which(name.endogenous==iVar)
                index.iVar <- which(rowSums(indexArr.lower.tri==iEndogenous)>0)

                ##  d sqrt(x) / d x = 1/(2 sqrt(x)) = sqrt(x) / (2*x)
                dOmega.dtheta[[iVar]][index.lower.tri[index.iVar]] <- var.coef["sigma2"]*dOmega.dtheta[["sigma2"]][index.lower.tri[index.iVar]]/(2*var.coef[iVar])
                dOmega.dtheta[[iVar]] <- symmetrize(dOmega.dtheta[[iVar]])
            }
            
            dimnames(dOmega.dtheta[[iVar]]) <- list(name.endogenous, name.endogenous)            
        }
    }
    
    ## *** correlation
    if("NULL" %in% class.cor == FALSE){
        for(iVar in name.corcoef){
            dOmega.dtheta[[iVar]] <- Msigma2.base0 * var.coef["sigma2"] * (M.corcoef==iVar)
        }
    }

### ** second order
    d2Omega.d2theta <- list()

    if(second.order){

        if("NULL" %in% class.var == FALSE){
            for(iVar in name.other){ ## iVar <- name.other[1]
                d2Omega.d2theta[["sigma"]][[iVar]] <- dOmega.dtheta[[iVar]]/var.coef["sigma2"]
            }
        }

        if("NULL" %in% class.cor == FALSE){
            for(iVar in name.corcoef){
                d2Omega.d2theta[["sigma"]][[iVar]] <- dOmega.dtheta[[iVar]]/var.coef["sigma2"]
            }
        }

        if("NULL" %in% class.var == FALSE && "NULL" %in% class.cor == FALSE){
            for(iVar1 in name.other){ ## iVar <- name.other[1]
                for(iVar2 in name.corcoef){
                    M.tempo <- matrix(1, nrow = n.endogenous, ncol = n.endogenous)
                    M.tempo[index.lower.tri] <- cor.coef[M.corcoef[index.lower.tri]]
                    M.tempo <- symmetrize(M.tempo, update.upper = TRUE)
                    
                    d2Omega.d2theta[[iVar1]][[iVar2]] <- dOmega.dtheta[[iVar1]]/M.tempo
                    diag(d2Omega.d2theta[[iVar1]][[iVar2]]) <- 0
                }
            }
        }

    }
    
### ** export
    return(list(dmu.dtheta = dmu.dtheta,
                dOmega.dtheta = dOmega.dtheta,
                d2Omega.d2theta = d2Omega.d2theta))
    
}

## * prepareScore2.lme
#' @rdname prepareScore2
#' @export
prepareScore2.lme <- function(object, X, 
                              param, attr.param,
                              second.order,
                              n.cluster, n.endogenous, name.endogenous, index.obs){

    resGLS <- prepareScore2.gls(object, X = X, param = param, attr.param = attr.param,
                                n.cluster = n.cluster,
                                second.order = second.order,
                                n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                index.obs = index.obs)
        
### ** random effect
    name.rancoef <- attr.param$ran.coef
    resGLS$dOmega.dtheta[[name.rancoef]] <- matrix(1, nrow = n.endogenous, ncol = n.endogenous,
                                                   dimnames = list(name.endogenous,name.endogenous))

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
                                 name.endogenous = NULL, name.latent = NULL,
                                 second.order = FALSE){

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
    
### ** Update first order partial derivatives with current estimates
    prepareScore2$dtheta <- skeletonDtheta(object, data = data,
                                           dt.param.all = prepareScore2$skeleton$dt.param,
                                           param2originalLink = prepareScore2$skeleton$param2originalLink,
                                           name.endogenous = name.endogenous, 
                                           name.latent = name.latent,
                                           B = prepareScore2$skeleton$value$B,
                                           alpha.XGamma = prepareScore2$skeleton$value$alpha.XGamma,
                                           Lambda = prepareScore2$skeleton$value$Lambda,
                                           Psi = prepareScore2$skeleton$value$Psi)

### ** Compute second order partial derivatives with current estimates
    if(second.order){
        prepareScore2$dtheta2 <- skeletonDtheta2(object, data = data,
                                                 dt.param.all = prepareScore2$skeleton$dt.param,
                                                 param2originalLink = prepareScore2$skeleton$param2originalLink,
                                                 name.endogenous = name.endogenous, 
                                                 name.latent = name.latent,
                                                 B = prepareScore2$skeleton$value$B,
                                                 alpha.XGamma = prepareScore2$skeleton$value$alpha.XGamma,
                                                 Lambda = prepareScore2$skeleton$value$Lambda,
                                                 Psi = prepareScore2$skeleton$value$Psi)
    }
    
### ** Export
    return(prepareScore2)    
    
}




#----------------------------------------------------------------------
### prepareScore2.R ends here
