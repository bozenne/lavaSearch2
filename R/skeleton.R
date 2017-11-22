### skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (10:35) 
## Version: 
## Last-Updated: nov 17 2017 (10:44) 
##           By: Brice Ozenne
##     Update #: 220
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - skeleton
#' @title
#' @description
#' @name skeleton
#' 
#' @param object a latent variable model
#' @param as.lava should the name of the links be used to name the parameters?
#' Otherwise uses the labels (when defined) of each parameter.
#' @param name.endogenous name of the endogenous variables
#' @param n.endogenous number of endogenous variables
#' @param name.latent name of the latent variables
#' @param n.latent number of latent variables
#' @param update.value should the design matrices include the current model estimates?
#' Otherwise put NA for each paramter to be estimated.
#' @param p [optional] vector of parameters at which to evaluate the score.
#' @param data [optional] data set.
#' @details
#' When the use specify names for the parameters (e.g. Y1[mu:sigma]) or uses constrains (Y1~beta*X1), \code{as.lava=FALSE} will use the names specified by the user (e.g. mu, sigma, beta) while \code{as.lava=TRUE} will use the name of the first link defining the parameter.
#'
#' @examples
#' ## without constrain
#' m <- lvm(Y1~X1+X2+eta,Y2~X3+eta,Y3~eta)
#' latent(m) <- ~eta
#' 
#' e <- estimate(m,sim(m,1e2))
#' skeleton(e, as.lava = TRUE,
#'          name.endogenous = endogenous(e), n.endogenous = 3,
#'          name.latent = latent(e), n.latent = 1,
#'          update.value = FALSE)
#' skeleton(e, as.lava = TRUE,
#'          name.endogenous = endogenous(e), n.endogenous = 3,
#'          name.latent = latent(e), n.latent = 1,
#'          update.value = TRUE,
#'          data = as.matrix(e$data$model.frame), p = pars(e))$value
#'
#' ## with constrains
#' m <- lvm(Y[mu:sigma] ~ beta*X1+X2)
#' e <- estimate(m, sim(m,1e2))
#'
#' skeleton(e, as.lava = TRUE,
#'          name.endogenous = "Y", n.endogenous = 1,
#'          name.latent = NULL, n.latent = 0,
#'          update.value = FALSE)$skeleton
#' 
#' skeleton(e, as.lava = FALSE,
#'          name.endogenous = "Y", n.endogenous = 1,
#'          name.latent = NULL, n.latent = 0,
#'          update.value = FALSE)$skeleton
#' 
#' @export
`skeleton` <-
    function(object, ...) UseMethod("skeleton")


## * skeleton.lvm
#' @rdname skeleton
#' @export
skeleton.lvm <- function(object, as.lava,
                         name.endogenous, name.latent){

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)

### ** prepare
    dt.param.all <- coefType(object, as.lava = FALSE)
    if(as.lava){
        param2originalLink <- dt.param.all[!is.na(lava)][,setNames(originalLink,param)]
    }else{
        param2originalLink <- dt.param.all[!is.na(param),setNames(param,param)]
    }
    
    skeleton <- list()
    value <- list()

### ** Measurement model
    ## *** nu
    dt.param.nu <-  dt.param.all[detail=="nu",.(value,param,Y),by="name"]
    value$nu <- setNames(dt.param.nu$value,dt.param.nu$Y)
    skeleton$nu <- setNames(param2originalLink[dt.param.nu$param],dt.param.nu$Y)

    ## *** X K
    dt.param.K <- dt.param.all[detail=="K",.(value,param,X),by="Y"]

    if(NROW(dt.param.K)>0){
        value$K <- setNames(lapply(1:n.endogenous, function(iEndogenous){
            dt.param.K$value[dt.param.K$Y==name.endogenous[iEndogenous]]
        }), name.endogenous)
            
        skeleton$K <- setNames(lapply(1:n.endogenous, function(iEndogenous){
            param2originalLink[dt.param.K$param[dt.param.K$Y==name.endogenous[iEndogenous]]]
        }), name.endogenous)
    
        skeleton$XK <- setNames(lapply(1:n.endogenous, function(iEndogenous){
            dt.param.K$X[dt.param.K$Y==name.endogenous[iEndogenous]]
        }), name.endogenous)
    }
    
    ## *** Lambda
    if(n.latent>0){
        ## define matrix
        value$Lambda <- matrix(0,nrow = n.latent, ncol = n.endogenous,
                               dimnames = list(name.latent,name.endogenous))
        skeleton$Lambda <- matrix(as.character(NA),nrow = n.latent, ncol = n.endogenous,
                                  dimnames = list(name.latent,name.endogenous))
        ## update according to the model
        dt.param.Lambda <- dt.param.all[detail %in% "Lambda",
                                       .(index = which(name.latent %in% X) + n.latent*(which(name.endogenous %in% Y)-1), param, value),
                                       by = name]
        
        skeleton$Lambda[dt.param.Lambda[is.na(value)]$index] <- setNames(param2originalLink[dt.param.Lambda[is.na(value)]$param],dt.param.Lambda$Y)
        value$Lambda[dt.param.Lambda[!is.na(value)]$index] <- setNames(dt.param.Lambda[!is.na(value)]$value,dt.param.Lambda$Y)
        value$Lambda[!is.na(skeleton$Lambda)] <- NA
    }
    
    ## *** Sigma    
    ## define matrix
    value$Sigma <- matrix(0,nrow = n.endogenous, ncol = n.endogenous,
                          dimnames = list(name.endogenous,name.endogenous))
    skeleton$Sigma <- matrix(as.character(NA),nrow = n.endogenous, ncol = n.endogenous,
                             dimnames = list(name.endogenous,name.endogenous))

    ## update according to the model
    dt.param.Sigma <- dt.param.all[detail %in% c("Sigma_var","Sigma_cov"),
                                   .(index = which(name.endogenous %in% X) + n.endogenous*(which(name.endogenous %in% Y)-1), param, value),
                                   by = name]
    
    skeleton$Sigma[dt.param.Sigma[is.na(value)]$index] <- param2originalLink[dt.param.Sigma[is.na(value)]$param]
    value$Sigma[dt.param.Sigma[!is.na(value)]$index] <- dt.param.Sigma[!is.na(value)]$value

    ## symmetrize
    skeleton$Sigma <- symmetrize(skeleton$Sigma, update.upper = TRUE)
    value$Sigma <- symmetrize(value$Sigma, update.upper = TRUE)
    value$Sigma[!is.na(skeleton$Sigma)] <- NA

### ** Structural model
    if(n.latent>0){
        ## *** alpha 
        dt.param.alpha <-  dt.param.all[detail=="alpha",.(value,param,Y),by="name"]
        value$alpha <- setNames(dt.param.alpha$value,dt.param.alpha$Y)
        skeleton$alpha <- param2originalLink[setNames(dt.param.alpha$param,dt.param.alpha$Y)]

        ## *** X Gamma
        dt.param.Gamma <- dt.param.all[detail=="Gamma",.(value,param,X),by="Y"]

        if(NROW(dt.param.Gamma)>0){
            value$Gamma <- setNames(lapply(1:n.latent, function(iLatent){
                dt.param.Gamma$value[dt.param.Gamma$Y==name.latent[iLatent]]
            }), name.latent)
            
            skeleton$Gamma <- setNames(lapply(1:n.latent, function(iLatent){
                param2originalLink[dt.param.Gamma$param[dt.param.Gamma$Y==name.latent[iLatent]]]
            }), name.latent)
    
            skeleton$XGamma <- setNames(lapply(1:n.latent, function(iLatent){
                dt.param.Gamma$X[dt.param.Gamma$Y==name.latent[iLatent]]
            }), name.latent)
        }
        
        ## *** B
        ## define matrix
        value$B <- matrix(0,nrow = n.latent, ncol = n.latent,
                          dimnames = list(name.latent,name.latent))
        skeleton$B <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                             dimnames = list(name.latent,name.latent))

        if(any("B" %in% dt.param.all$detail)){
            ## update according to the model
            dt.param.B <- dt.param.all[detail == "B", .(index = which(name.latent %in% X) + n.latent*(which(name.latent %in% Y)-1), param, value), by = name]
            skeleton$B[dt.param.B[is.na(value)]$index] <- param2originalLink[dt.param.B[is.na(value),param]]
            value$B[dt.param.B[!is.na(value)]$index] <- dt.param.B[!is.na(value),value]
            value$B[!is.na(skeleton$B)] <- NA            
        }
    
        ## *** Psi    
        ## define matrix
        value$Psi <- matrix(0,nrow = n.latent, ncol = n.latent,
                            dimnames = list(name.latent,name.latent))
        skeleton$Psi <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                               dimnames = list(name.latent,name.latent))
    
        ## update according to the model
        dt.param.Psi <- dt.param.all[detail %in% c("Psi_var","Psi_cov"), .(index = which(name.latent %in% X) + n.latent*(which(name.latent %in% Y)-1), param, value, Y), by = name]
        skeleton$Psi[dt.param.Psi[is.na(value)]$index] <- param2originalLink[dt.param.Psi[is.na(value)]$param]
        value$Psi[dt.param.Psi[!is.na(value)]$index] <- dt.param.Psi[!is.na(value)]$value

        ## symmetrize
        skeleton$Psi <- symmetrize(skeleton$Psi, update.upper = TRUE)
        value$Psi <- symmetrize(value$Psi, update.upper = TRUE)
        value$Psi[!is.na(skeleton$Psi)] <- NA
    }
    
### ** export
    return(list(skeleton = skeleton,
                value = value,
                dt.param = dt.param.all,
                param2originalLink = param2originalLink)
           )
}


## * skeleton.lvmfit
#' @rdname skeleton
#' @export
skeleton.lvmfit <- function(object, as.lava, 
                            p, data,
                            name.endogenous, name.latent){

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
    n <- NROW(data)

### ** Compute skeleton
    if(is.null(object$prepareScore2$skeleton)){
        OS <- skeleton(lava::Model(object), as.lava = as.lava,
                       name.endogenous = name.endogenous, name.latent = name.latent)
    }else{
        OS <- object$prepareScore2$skeleton
    }

### ** Update skeleton with the current values
    
    ## *** nu
    index.update <- which(!is.na(OS$skeleton$nu))
    OS$value$nu[index.update] <- p[OS$skeleton$nu[index.update]]

    ## *** K 
    for(iY in 1:n.endogenous){ # iY <- 3
        if(length(OS$skeleton$K[[iY]])>0){
            index.update <- which(!is.na(OS$skeleton$K[[iY]]))
            OS$value$K[[iY]][index.update] <- p[OS$skeleton$K[[iY]][index.update]]
        }
    }

    ## *** Lambda
    if(n.latent>0){
        index.update <- which(!is.na(OS$skeleton$Lambda))
        OS$value$Lambda[index.update] <- p[OS$skeleton$Lambda[index.update]]
    }
    
    ## *** Sigma
    index.update <- which(!is.na(OS$skeleton$Sigma))
    OS$value$Sigma[index.update] <- p[OS$skeleton$Sigma[index.update]]

    ## *** linear predictor
    OS$value$nu.XK <- matrix(NA,nrow = n, ncol = n.endogenous, byrow = TRUE,
                             dimnames = list(NULL,name.endogenous))
    for(iY in 1:n.endogenous){ # iY <- 1
        iY2 <- name.endogenous[iY]
        if(length(OS$value$K[[iY2]])>0){
            OS$value$nu.XK[,iY2] <- OS$value$nu[iY2] + data[,OS$skeleton$XK[[iY2]],drop=FALSE] %*% OS$value$K[[iY2]]
        }else{
            OS$value$nu.XK[,iY2] <- OS$value$nu[iY2]
        }
    }
        
### ** Structural model
    if(n.latent>0){
        ## *** alpha
        index.update <- which(!is.na(OS$skeleton$alpha))
        OS$value$alpha[index.update] <- p[OS$skeleton$alpha[index.update]]

        ## *** Gamma
        for(iLatent in 1:n.latent){
            if(length(OS$skeleton$Gamma[[iLatent]])>0){
                index.update <- which(!is.na(OS$skeleton$Gamma[[iLatent]]))
                OS$value$Gamma[[iLatent]][index.update] <- p[OS$skeleton$Gamma[[iLatent]][index.update]]
            }
        }
        
        ## *** B
        index.update <- which(!is.na(OS$skeleton$B))
        OS$value$B[index.update] <- p[OS$skeleton$B[index.update]]

        ## *** Psi
        index.update <- which(!is.na(OS$skeleton$Psi))
        OS$value$Psi[index.update] <- p[OS$skeleton$Psi[index.update]]

        ## *** linear predictor
        OS$value$alpha.XGamma <- matrix(NA,nrow = n, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL,name.latent))
        for(iLatent in 1:n.latent){
            iLatent2 <- name.latent[iLatent]
                if(length(OS$value$Gamma[[iLatent2]])>0){
                    OS$value$alpha.XGamma[,iLatent2] <- OS$value$alpha[iLatent2] + data[,OS$skeleton$XGamma[[iLatent2]],drop=FALSE] %*% OS$value$Gamma[[iLatent2]]
                }else{
                    OS$value$alpha.XGamma[,iLatent2] <- OS$value$alpha[iLatent2]
                }
        }
    }
           
### ** Export
    return(OS)
}


## * skeletonDtheta
#' @rdname skeleton
#' @export
`skeletonDtheta` <-
    function(object, ...) UseMethod("skeletonDtheta")

## * skeletonDtheta.lvm
#' @rdname skeleton
#' @export
skeletonDtheta.lvm <- function(object, data,
                               dt.param.all, param2originalLink,
                               name.endogenous, name.latent, ...){

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
    
    dt.param <- dt.param.all[is.na(value) & marginal == FALSE & factice == FALSE]
    Utype.by.detail <- dt.param[,.(n.type = length(unique(detail))), by = param][["n.type"]]
    if(any(Utype.by.detail>1)){
        stop("cannot constrain two parameters of different types to be equal \n")
    }
    name.param <- dt.param[!duplicated(param),param]
    n.param <- length(name.param)

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
        type = type,
        toUpdate = toUpdate
    ))
}


## * skeletonDtheta.lvmfit
#' @rdname skeleton
#' @export
skeletonDtheta.lvmfit <- function(object, data,
                                  dt.param.all, param2originalLink,
                                  name.endogenous, name.latent,
                                  B, alpha.XGamma, Lambda, Psi,
                                  ...){

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
  
### ** Initialize partial derivatives
    if(is.null(object$prepareScore2$dtheta)){
        OD <- skeletonDtheta(lava::Model(object), data = data,
                             dt.param.all = dt.param.all,
                             param2originalLink = param2originalLink,
                             name.endogenous = name.endogenous, 
                             name.latent = name.latent)
    }else{
        OD <- object$prepareScore2$dtheta
    }

### ** Update partial derivatives
    
    
    if(any(OD$toUpdate)){
        type2update <- OD$type[OD$toUpdate]
            
        OD$iIB <- solve(diag(1,n.latent,n.latent)-B)
        OD$alpha.XGamma.iIB <- alpha.XGamma %*% OD$iIB
        OD$iIB.Lambda <-  OD$iIB %*% Lambda    
        OD$Psi.iIB <- Psi %*% OD$iIB
        OD$tLambda.tiIB.Psi.iIB <- t(OD$iIB.Lambda) %*% OD$Psi.iIB
        
        ## *** mean parameters
        type.meanparam <- type2update[type2update %in% c("alpha","Lambda","Gamma","B")]
        n.meanparam <- length(type.meanparam)
        name.meanparam <- names(type.meanparam)

        if(n.meanparam>0){
            for(iP in 1:n.meanparam){ # iP <- 1
                iType <- type.meanparam[iP]
                iName <- name.meanparam[iP]
            
                if(iType == "alpha"){
                    OD$dmu.dtheta[[iName]] <- OD$dmu.dtheta[[iName]] %*% OD$iIB.Lambda
                }else if(iType == "Gamma"){
                    OD$dmu.dtheta[[iName]] <- OD$dmu.dtheta[[iName]] %*% OD$iIB.Lambda 
                }else if(iType == "Lambda"){
                    OD$dmu.dtheta[[iName]] <- OD$alpha.XGamma.iIB %*% OD$dLambda.dtheta[[iName]]
                }else if(iType == "B"){
                    OD$dmu.dtheta[[iName]] <- OD$alpha.XGamma.iIB %*% OD$dB.dtheta[[iName]] %*% OD$iIB.Lambda
                }

                colnames(OD$dmu.dtheta[[iName]]) <- name.endogenous
            }
        }

        ## *** variance-covariance parameters
        type.vcovparam <- type2update[type2update %in% c("Psi_var","Psi_cov","Lambda","B")]
        n.vcovparam <- length(type.vcovparam)
        name.vcovparam <- names(type.vcovparam)

        if(n.vcovparam>0){
            for(iP in 1:n.vcovparam){ # iP <- 1
                iType <- type.vcovparam[iP]
                iName <- name.vcovparam[iP]
        
                if(iType %in% "Psi_var"){
                    OD$dOmega.dtheta[[iName]] <-  t(OD$iIB.Lambda) %*% OD$dPsi.dtheta[[iName]] %*% OD$iIB.Lambda
                }else if(iType %in% "Psi_cov"){
                    OD$dOmega.dtheta[[iName]] <-  t(OD$iIB.Lambda) %*% OD$dPsi.dtheta[[iName]] %*% OD$iIB.Lambda
                }else if(iType == "Lambda"){
                    OD$dOmega.dtheta[[iName]] <- OD$tLambda.tiIB.Psi.iIB %*% OD$dLambda.dtheta[[iName]]
                    OD$dOmega.dtheta[[iName]] <- OD$dOmega.dtheta[[iName]] + t(OD$dOmega.dtheta[[iName]])
                }else if(iType == "B"){
                    OD$dOmega.dtheta[[iName]] <- OD$tLambda.tiIB.Psi.iIB %*% OD$dB.dtheta[[iName]] %*% OD$iIB.Lambda
                    OD$dOmega.dtheta[[iName]] <- OD$dOmega.dtheta[[iName]] + t(OD$dOmega.dtheta[[iName]])
                }

                colnames(OD$dOmega.dtheta[[iName]]) <- name.endogenous
                rownames(OD$dOmega.dtheta[[iName]]) <- name.endogenous
            }
        }
        
    }

### ** Export
    return(OD)

}

## * skeletonDtheta2
#' @rdname skeleton
#' @export
`skeletonDtheta2` <-
    function(object, ...) UseMethod("skeletonDtheta2")

## * skeletonDtheta2.lvmfit
#' @rdname skeleton
#' @export
skeletonDtheta2.lvmfit <- function(object, data,
                                   dt.param.all, param2originalLink,
                                   name.endogenous, name.latent,
                                   B, alpha.XGamma, Lambda, Psi,
                                   ...){

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
  
### ** second order partial derivatives
    browser()    
    if(any(OD$toUpdate)){
        type2update <- OD$type[OD$toUpdate]
            
        OD$iIB <- solve(diag(1,n.latent,n.latent)-B)
        OD$alpha.XGamma.iIB <- alpha.XGamma %*% OD$iIB
        OD$iIB.Lambda <-  OD$iIB %*% Lambda    
        OD$Psi.iIB <- Psi %*% OD$iIB
        OD$tLambda.tiIB.Psi.iIB <- t(OD$iIB.Lambda) %*% OD$Psi.iIB
        
        ## *** mean parameters
        type.meanparam <- type2update[type2update %in% c("alpha","Lambda","Gamma","B")]
        n.meanparam <- length(type.meanparam)
        name.meanparam <- names(type.meanparam)

        if(n.meanparam>0){
            for(iP in 1:n.meanparam){ # iP <- 1
                iType <- type.meanparam[iP]
                iName <- name.meanparam[iP]
            
                if(iType == "alpha"){
                    OD$dmu.dtheta[[iName]] <- OD$dmu.dtheta[[iName]] %*% OD$iIB.Lambda
                }else if(iType == "Gamma"){
                    OD$dmu.dtheta[[iName]] <- OD$dmu.dtheta[[iName]] %*% OD$iIB.Lambda 
                }else if(iType == "Lambda"){
                    OD$dmu.dtheta[[iName]] <- OD$alpha.XGamma.iIB %*% OD$dLambda.dtheta[[iName]]
                }else if(iType == "B"){
                    OD$dmu.dtheta[[iName]] <- OD$alpha.XGamma.iIB %*% OD$dB.dtheta[[iName]] %*% OD$iIB.Lambda
                }

                colnames(OD$dmu.dtheta[[iName]]) <- name.endogenous
            }
        }

        ## *** variance-covariance parameters
        type.vcovparam <- type2update[type2update %in% c("Psi_var","Psi_cov","Lambda","B")]
        n.vcovparam <- length(type.vcovparam)
        name.vcovparam <- names(type.vcovparam)

        if(n.vcovparam>0){
            for(iP in 1:n.vcovparam){ # iP <- 1
                iType <- type.vcovparam[iP]
                iName <- name.vcovparam[iP]
        
                if(iType %in% "Psi_var"){
                    OD$dOmega.dtheta[[iName]] <-  t(OD$iIB.Lambda) %*% OD$dPsi.dtheta[[iName]] %*% OD$iIB.Lambda
                }else if(iType %in% "Psi_cov"){
                    OD$dOmega.dtheta[[iName]] <-  t(OD$iIB.Lambda) %*% OD$dPsi.dtheta[[iName]] %*% OD$iIB.Lambda
                }else if(iType == "Lambda"){
                    OD$dOmega.dtheta[[iName]] <- OD$tLambda.tiIB.Psi.iIB %*% OD$dLambda.dtheta[[iName]]
                    OD$dOmega.dtheta[[iName]] <- OD$dOmega.dtheta[[iName]] + t(OD$dOmega.dtheta[[iName]])
                }else if(iType == "B"){
                    OD$dOmega.dtheta[[iName]] <- OD$tLambda.tiIB.Psi.iIB %*% OD$dB.dtheta[[iName]] %*% OD$iIB.Lambda
                    OD$dOmega.dtheta[[iName]] <- OD$dOmega.dtheta[[iName]] + t(OD$dOmega.dtheta[[iName]])
                }

                colnames(OD$dOmega.dtheta[[iName]]) <- name.endogenous
                rownames(OD$dOmega.dtheta[[iName]]) <- name.endogenous
            }
        }
        
    }

### ** Export
    return(OD)

}


##----------------------------------------------------------------------
### skeleton.R ends here

