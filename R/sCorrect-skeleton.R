### sCorrect-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (10:35) 
## Version: 
## Last-Updated: jan  7 2020 (14:35) 
##           By: Brice Ozenne
##     Update #: 1474
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - skeleton
#' @title Pre-computation for the Score
#' @description Pre-compute quantities that are necessary to compute the score of a lvm model.
#' @name skeleton
#' 
#' @param object a \code{lvm} object.
#' @param X [matrix] design matrix containing the covariates for each endogeneous and latent variable.
#' @param endogeneous [character vector] the name of the endogeneous variables.
#' @param latent [character vector] the name of the latent variables.
#' @param ... [internal] only used by the generic method.
#' 
#' @details
#' When the use specify names for the coefficients (e.g. Y1[mu:sigma]) or uses constrains (Y1~beta*X1), \code{as.lava=FALSE} will use the names specified by the user (e.g. mu, sigma, beta) while \code{as.lava=TRUE} will use the name of the first link defining the coefficient.
#'
#' @examples
#' \dontrun{
#' skeleton <- lavaSearch2::skeleton
#' skeleton.lvm <- lavaSearch2::skeleton.lvm
#' skeleton.lvmfit <- lavaSearch2::skeleton.lvmfit
#' 
#' ## without constrain
#' m <- lvm(Y1~X1+X2+eta,Y2~X3+eta,Y3~eta)
#' latent(m) <- ~eta
#' 
#' e <- estimate(m, lava::sim(m,1e2))
#' M.data <- as.matrix(model.frame(e))
#'
#' skeleton(e$model, as.lava = TRUE,
#'          name.endogenous = endogenous(e), n.endogenous = 3,
#'          name.latent = latent(e), 
#'          update.value = FALSE)
#' skeleton(e, data = M.data, p = pars(e), as.lava = TRUE,
#'          name.endogenous = endogenous(e), n.endogenous = 3,
#'          name.latent = latent(e), 
#'          update.value = TRUE)
#'
#' ## with constrains
#' m <- lvm(Y[mu:sigma] ~ beta*X1+X2)
#' e <- estimate(m, lava::sim(m,1e2))
#' M.data <- as.matrix(model.frame(e))
#'
#' skeleton(e$model, as.lava = TRUE,
#'          name.endogenous = "Y", n.endogenous = 1,
#'          name.latent = NULL, 
#'          update.value = FALSE)$skeleton
#' 
#' skeleton(e, data = M.data, p = pars(e), as.lava = FALSE,
#'          name.endogenous = "Y", n.endogenous = 1,
#'          name.latent = NULL, 
#'          update.value = FALSE)$skeleton
#' 
#'}
#' @concept small sample inference
#' @concept derivative of the score equation
#' @keywords internal


## * skeleton
#' @rdname skeleton
#' @export
skeleton <- function(object, X,
                     endogenous, latent,
                     n.cluster, index.Omega){
    if(lava.options()$debug){cat("skeleton \n")}
    
    n.endogenous <- length(endogenous)
    n.latent <- length(latent)
    n.obs <- NROW(X)
    obsByEndoInX <- tapply(1:n.obs,X$XXendogenousXX,list)

    ## ** extract table type
    type <- coefType(object, as.lava = FALSE)
    type.theta <- type[type$marginal == FALSE,,drop=FALSE]
    theta.value <- list()
    theta.param <- list()

    ## ** prepare value and param
    if("nu" %in% type.theta$detail){
        type.nu <- type.theta[type.theta$detail %in% "nu",,drop=FALSE]

        theta.value$nu <- setNames(rep(NA, n.endogenous), endogenous)
        if(any(!is.na(type.nu$value))){
            theta.value$nu[type.nu[!is.na(type.nu$value),"Y"]] <- type.nu[!is.na(type.nu$value),"value"]
        }

        theta.param$nu <- setNames(rep(as.character(NA), n.endogenous), endogenous)
        if(any(!is.na(type.nu$param))){
            theta.param$nu[type.nu[!is.na(type.nu$param),"Y"]] <- type.nu[!is.na(type.nu$param),"param"]
        }

        theta.param$Xnu <- matrix(NA,
                                  nrow = n.cluster, ncol = n.endogenous, byrow = TRUE,
                                  dimnames = list(NULL,endogenous))
        for(iC in 1:n.cluster){
            theta.param$Xnu[iC,index.Omega[[iC]]] <- 1
        }        
    }

    if("alpha" %in% type.theta$detail){
        type.alpha <- type.theta[type.theta$detail %in% "alpha",,drop=FALSE]
        theta.value$alpha <- setNames(rep(NA, n.latent), latent)
        theta.param$alpha <- setNames(rep(as.character(NA), n.latent), latent)
        theta.param$Xalpha <- rep(1, times = n.cluster)
        
        if(any(!is.na(type.alpha$value))){
            theta.value$alpha[type.alpha[!is.na(type.alpha$value),"Y"]] <- type.alpha[!is.na(type.alpha$value),"value"]
        }
        if(any(!is.na(type.alpha$param))){
            theta.param$alpha[type.alpha[!is.na(type.alpha$param),"Y"]] <- type.alpha[!is.na(type.alpha$param),"param"]
        }
    }

    if("K" %in% type.theta$detail){
        type.K <- type.theta[type.theta$detail %in% "K",,drop=FALSE]
        K.exogenous <- unique(type.K$X)
        
        theta.value$K <- matrix(0, nrow = length(K.exogenous), ncol = n.endogenous,
                                dimnames = list(K.exogenous, endogenous))
        theta.param$K <- matrix(as.character(NA), nrow = length(K.exogenous), ncol = n.endogenous,
                                dimnames = list(K.exogenous, endogenous))
        theta.param$XK <- lapply(endogenous, function(iE){ ## iE <- endogenous[1]
            iXK <- matrix(NA, nrow = n.cluster, ncol = length(K.exogenous),
                          dimnames = list(NULL, K.exogenous))
            iXK[X[obsByEndoInX[[iE]],"XXclusterXX"],] <- as.matrix(X[obsByEndoInX[[iE]],K.exogenous,drop=FALSE])
            return(iXK)
        })
        names(theta.param$XK) <- endogenous
        
        for(iK in 1:NROW(type.K)){ ## iK <- 1
            theta.value$K[type.K[iK,"X"],type.K[iK,"Y"]] <- type.K[iK,"value"]
            theta.param$K[type.K[iK,"X"],type.K[iK,"Y"]] <- type.K[iK,"param"]
        }
    }

    if("Gamma" %in% type.theta$detail){
        type.Gamma <- type.theta[type.theta$detail %in% "Gamma",,drop=FALSE]
        Gamma.exogenous <- unique(type.Gamma$X)
        
        theta.value$Gamma <- matrix(0, nrow = length(Gamma.exogenous), ncol = n.latent,
                                    dimnames = list(Gamma.exogenous, latent))
        theta.param$Gamma <- matrix(as.character(NA), nrow = length(Gamma.exogenous), ncol = n.latent,
                                    dimnames = list(Gamma.exogenous, latent))
        theta.param$XGamma <- lapply(obsByEndoInX[latent], function(iIndex){as.matrix(X[iIndex,Gamma.exogenous,drop=FALSE])})

        for(iGamma in 1:NROW(type.Gamma)){ ## iGamma <- 1
            theta.value$Gamma[type.Gamma[iGamma,"X"],type.Gamma[iGamma,"Y"]] <- type.Gamma[iGamma,"value"]
            theta.param$Gamma[type.Gamma[iGamma,"X"],type.Gamma[iGamma,"Y"]] <- type.Gamma[iGamma,"param"]
        }
    }

    if("Lambda" %in% type.theta$detail){
        type.Lambda <- type.theta[type.theta$detail %in% "Lambda",,drop=FALSE]
        
        theta.value$Lambda <- matrix(0, nrow = n.latent, ncol = n.endogenous,
                                     dimnames = list(latent, endogenous))
        theta.param$Lambda <- matrix(as.character(NA), nrow = n.latent, ncol = n.endogenous,
                                     dimnames = list(latent, endogenous))

        for(iLambda in 1:NROW(type.Lambda)){ ## iLambda <- 1
            theta.value$Lambda[type.Lambda[iLambda,"X"],type.Lambda[iLambda,"Y"]] <- type.Lambda[iLambda,"value"]
            theta.param$Lambda[type.Lambda[iLambda,"X"],type.Lambda[iLambda,"Y"]] <- type.Lambda[iLambda,"param"]
        }
    }

    if("B" %in% type.theta$detail){
        type.B <- type.theta[type.theta$detail %in% "B",,drop=FALSE]
        
        theta.value$B <- matrix(0, nrow = n.latent, ncol = n.latent,
                                     dimnames = list(latent, latent))
        theta.param$B <- matrix(as.character(NA), nrow = n.latent, ncol = n.latent,
                                     dimnames = list(latent, latent))

        for(iB in 1:NROW(type.B)){ ## iB <- 1
            theta.value$B[type.B[iB,"X"],type.B[iB,"Y"]] <- type.B[iB,"value"]
            theta.param$B[type.B[iB,"X"],type.B[iB,"Y"]] <- type.B[iB,"param"]
        }
    }

    if(any(c("Sigma_var", "Sigma_cov") %in% type.theta$detail)){
        type.Sigma_var <- type.theta[type.theta$detail %in% "Sigma_var",,drop=FALSE]
        type.Sigma_cov <- type.theta[type.theta$detail %in% "Sigma_cov",,drop=FALSE]

        theta.value$Sigma <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                                     dimnames = list(endogenous, endogenous))
        theta.param$Sigma <- matrix(as.character(NA), nrow = n.endogenous, ncol = n.endogenous,
                                     dimnames = list(endogenous, endogenous))

        if(NROW(type.Sigma_var)>0){
            for(iSigma_var in 1:NROW(type.Sigma_var)){ ## iSigma_var <- 1
                theta.value$Sigma[type.Sigma_var[iSigma_var,"X"],type.Sigma_var[iSigma_var,"Y"]] <- type.Sigma_var[iSigma_var,"value"]
                theta.param$Sigma[type.Sigma_var[iSigma_var,"X"],type.Sigma_var[iSigma_var,"Y"]] <- type.Sigma_var[iSigma_var,"param"]
            }
        }
        if(NROW(type.Sigma_cov)>0){
            for(iSigma_cov in 1:NROW(type.Sigma_cov)){ ## iSigma_cov <- 1
                theta.value$Sigma[type.Sigma_cov[iSigma_cov,"X"],type.Sigma_cov[iSigma_cov,"Y"]] <- type.Sigma_cov[iSigma_cov,"value"]
                theta.value$Sigma[type.Sigma_cov[iSigma_cov,"Y"],type.Sigma_cov[iSigma_cov,"X"]] <- type.Sigma_cov[iSigma_cov,"value"]
                theta.param$Sigma[type.Sigma_cov[iSigma_cov,"X"],type.Sigma_cov[iSigma_cov,"Y"]] <- type.Sigma_cov[iSigma_cov,"param"]
                theta.param$Sigma[type.Sigma_cov[iSigma_cov,"Y"],type.Sigma_cov[iSigma_cov,"X"]] <- type.Sigma_cov[iSigma_cov,"param"]
            }
        }
    }else if(any(c("sigma2", "sigma2k", "cor") %in% type.theta$detail)){

        name.layer <- NULL
        if("sigma2" %in% type.theta$detail){
            name.layer <- c(name.layer,"sigma2")
        }
        if("sigma2k" %in% type.theta$detail){
            name.layer <- c(name.layer,"sigma2kX","sigma2kY")
        }
        if("cor" %in% type.theta$detail){
            name.layer <- c(name.layer,"cor")
        }
        n.layer <- length(name.layer)
        theta.param$SigmaValue <- array(0, dim = c(n.endogenous, n.endogenous, n.layer),
                                   dimnames = list(endogenous, endogenous, name.layer))
        theta.param$SigmaParam <- array(as.character(NA), dim = c(n.endogenous, n.endogenous, n.layer),
                                   dimnames = list(endogenous, endogenous, name.layer))

        if("sigma2" %in% type.theta$detail){
            type.sigma2 <- type.theta[type.theta$detail %in% "sigma2",,drop=FALSE]
            type.sigmacor2 <- type.theta[type.theta$detail %in% c("cor","sigma2"),,drop=FALSE]

            for(isigma2 in 1:NROW(type.sigmacor2)){ ## isigma2 <- 1
                theta.param$SigmaValue[type.sigmacor2[isigma2,"X"],type.sigmacor2[isigma2,"Y"],"sigma2"] <- type.sigma2[1,"value"]
                theta.param$SigmaValue[type.sigmacor2[isigma2,"Y"],type.sigmacor2[isigma2,"X"],"sigma2"] <- type.sigma2[1,"value"]
                
                theta.param$SigmaParam[type.sigmacor2[isigma2,"X"],type.sigmacor2[isigma2,"Y"],"sigma2"] <- type.sigma2[1,"param"]
                theta.param$SigmaParam[type.sigmacor2[isigma2,"Y"],type.sigmacor2[isigma2,"X"],"sigma2"] <- type.sigma2[1,"param"]
            }
        }
            

        if("sigma2k" %in% type.theta$detail){
            type.sigma2k <- type.theta[type.theta$detail %in% "sigma2k",,drop=FALSE]
            type.sigmacor2k <- type.theta[type.theta$detail %in% c("sigma2k","cor"),,drop=FALSE]
            
            for(isigma2k in 1:NROW(type.sigmacor2k)){ ## isigma2k <- 1
                theta.param$SigmaValue[type.sigmacor2k[isigma2k,"X"],type.sigmacor2k[isigma2k,"Y"],"sigma2kX"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"X"],"value"]
                theta.param$SigmaValue[type.sigmacor2k[isigma2k,"Y"],type.sigmacor2k[isigma2k,"X"],"sigma2kX"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"X"],"value"]

                theta.param$SigmaValue[type.sigmacor2k[isigma2k,"X"],type.sigmacor2k[isigma2k,"Y"],"sigma2kY"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"Y"],"value"]
                theta.param$SigmaValue[type.sigmacor2k[isigma2k,"Y"],type.sigmacor2k[isigma2k,"X"],"sigma2kY"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"Y"],"value"]

                theta.param$SigmaParam[type.sigmacor2k[isigma2k,"X"],type.sigmacor2k[isigma2k,"Y"],"sigma2kX"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"X"],"param"]
                theta.param$SigmaParam[type.sigmacor2k[isigma2k,"Y"],type.sigmacor2k[isigma2k,"X"],"sigma2kX"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"X"],"param"]

                theta.param$SigmaParam[type.sigmacor2k[isigma2k,"X"],type.sigmacor2k[isigma2k,"Y"],"sigma2kY"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"Y"],"param"]
                theta.param$SigmaParam[type.sigmacor2k[isigma2k,"Y"],type.sigmacor2k[isigma2k,"X"],"sigma2kY"] <- type.sigma2k[type.sigma2k$Y == type.sigmacor2k[isigma2k,"Y"],"param"]

            }
        }

        if("cor" %in% type.theta$detail){
            type.cor <- type.theta[type.theta$detail %in% "cor",,drop=FALSE]
            diag(theta.param$SigmaValue[,,"cor"]) <- 1
            for(icor in 1:NROW(type.cor)){ ## icor <- 1
                theta.param$SigmaValue[type.cor[icor,"X"],type.cor[icor,"Y"],"cor"] <- type.cor[icor,"value"]
                theta.param$SigmaValue[type.cor[icor,"Y"],type.cor[icor,"X"],"cor"] <- type.cor[icor,"value"]
                
                theta.param$SigmaParam[type.cor[icor,"X"],type.cor[icor,"Y"],"cor"] <- type.cor[icor,"param"]
                theta.param$SigmaParam[type.cor[icor,"Y"],type.cor[icor,"X"],"cor"] <- type.cor[icor,"param"]
            }
        }
    }    


    if(any(c("Psi_var", "Psi_cov") %in% type.theta$detail)){
        type.Psi_var <- type.theta[type.theta$detail %in% "Psi_var",,drop=FALSE]
        type.Psi_cov <- type.theta[type.theta$detail %in% "Psi_cov",,drop=FALSE]
        
        theta.value$Psi <- matrix(0, nrow = n.latent, ncol = n.latent,
                                     dimnames = list(latent, latent))
        theta.param$Psi <- matrix(as.character(NA), nrow = n.latent, ncol = n.latent,
                                     dimnames = list(latent, latent))

        if(NROW(type.Psi_var)>0){
            for(iPsi_var in 1:NROW(type.Psi_var)){ ## iPsi_var <- 1
                theta.value$Psi[type.Psi_var[iPsi_var,"X"],type.Psi_var[iPsi_var,"Y"]] <- type.Psi_var[iPsi_var,"value"]
                theta.param$Psi[type.Psi_var[iPsi_var,"X"],type.Psi_var[iPsi_var,"Y"]] <- type.Psi_var[iPsi_var,"param"]
            }
        }
        if(NROW(type.Psi_cov)>0){
            for(iPsi_cov in 1:NROW(type.Psi_cov)){ ## iPsi_cov <- 1
                theta.value$Psi[type.Psi_cov[iPsi_cov,"X"],type.Psi_cov[iPsi_cov,"Y"]] <- type.Psi_cov[iPsi_cov,"value"]
                theta.value$Psi[type.Psi_cov[iPsi_cov,"Y"],type.Psi_cov[iPsi_cov,"X"]] <- type.Psi_cov[iPsi_cov,"value"]
                theta.param$Psi[type.Psi_cov[iPsi_cov,"X"],type.Psi_cov[iPsi_cov,"Y"]] <- type.Psi_cov[iPsi_cov,"param"]
                theta.param$Psi[type.Psi_cov[iPsi_cov,"Y"],type.Psi_cov[iPsi_cov,"X"]] <- type.Psi_cov[iPsi_cov,"param"]
            }
        }
    }

    ## ** original link
    type.originalLink <- type[!is.na(type$originalLink),,drop=FALSE]
    originalLink2param <- setNames(type.originalLink$param,type.originalLink$originalLink)
    if(inherits(object,"lvm")){
        originalLink2param <- originalLink2param[coef(object)]
    }else if(inherits(object,"lvmfit")){
        originalLink2param <- originalLink2param[names(coef(object))]
    }
    
    ## ** type of parameters
    type.param <- type[!is.na(type$param),,drop=FALSE]
    type.mean <- c("nu","alpha","K","Gamma","Lambda","B")
    type.var <- c("Lambda","B","Sigma_var","Sigma_cov","sigma2","sigma2k","cor","Psi_var","Psi_cov")

    Uparam <- as.character(originalLink2param)
    Uparam.mean <- unique(type.param[type.param$detail %in% type.mean,"param"])
    Uparam.variance <- unique(type.param[type.param$detail %in% type.var,"param"])
    
    ## ** to update
    toUpdate <- c("nu" = "nu" %in% type.param$detail,
                  "K" = "K" %in% type.param$detail,
                  "Lambda" = "Lambda" %in% type.param$detail,
                  "Sigma" = ("Sigma_cov" %in% type.param$detail) || ("Sigma_var" %in% type.param$detail),
                  "SigmaValue" = ("sigma2" %in% type.param$detail) || ("sigma2k" %in% type.param$detail) || ("cor" %in% type.param$detail), 
                  "alpha" = "alpha" %in% type.param$detail,
                  "Gamma" = "Gamma" %in% type.param$detail,
                  "B" = "B" %in% type.param$detail,
                  "Psi" = ("Psi_cov" %in% type.param$detail) || ("Psi_var" %in% type.param$detail)
                  )    

    ## ** residuals
    theta.param$endogenous <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                                    dimnames = list(NULL, endogenous))
    for(iEndo in 1:length(endogenous)){ ## iEndo <- 1
        iIndexLong <- obsByEndoInX[[endogenous[iEndo]]]
        iIndexWide <- X[iIndexLong,"XXclusterXX"]
        theta.param$endogenous[iIndexWide,endogenous[iEndo]] <- X[iIndexLong,"XXvalueXX"]
    }

    ## ** export
    return(list(param = theta.param,
                value = theta.value,
                type = type,
                Uparam = Uparam,
                Uparam.mean = Uparam.mean,
                Uparam.variance = Uparam.variance,
                toUpdate.moment = toUpdate,
                originalLink2param = originalLink2param,
                obsByEndoInX = obsByEndoInX)
           )
}

## * skeletonDtheta
#' @rdname skeleton
skeletonDtheta <- function(object, X,
                           endogenous, latent,
                           missing.pattern, unique.pattern, name.pattern,
                           n.cluster, index.Omega){
    if(lava.options()$debug){cat("skeletonDtheta \n")}

    n.endogenous <- length(endogenous)
    n.latent <- length(latent)
    type <- object$type

    type <- type[!is.na(type$param),]
    name.param <- unique(type$param)
    n.param <- length(name.param)

    ## ** Compute partial derivative with respect to the matrices of parameters
    dOmega.dparam <- list() ## derivative regarding the variance
    dmat.dparam <- list() ## derivative regarding the matrix of parameters

    if("nu" %in% type$detail){
        type.nu <- type[type$detail == "nu",]
        Utype.nu <- unique(type.nu$param)
        nUtype.nu <- length(Utype.nu)

        dmat.dparam$nu <- setNames(vector(mode = "list", length = nUtype.nu), Utype.nu)
        
        for(iNu in Utype.nu){ ## iNu <- Utype.nu[1]
            dnu.dparam <- as.numeric(endogenous %in% type.nu[type.nu$param==iNu,"Y"])

            dmat.dparam$nu[[iNu]] <- matrix(NA,
                                            nrow = n.cluster, ncol = n.endogenous, byrow = TRUE,
                                            dimnames = list(NULL,endogenous))
            for(iP in name.pattern){ ## iP <- name.pattern[1]
                iIndex <- missing.pattern[[iP]]
                iY <- which(unique.pattern[iP,]==1)
                dmat.dparam$nu[[iNu]][iIndex,iY] <- rep(1, times = length(iIndex)) %o% dnu.dparam[iY]
            }
        }
    }

    if("alpha" %in% type$detail){
        type.alpha <- type[type$detail == "alpha",]
        Utype.alpha <- unique(type.alpha$param)
        nUtype.alpha <- length(Utype.alpha)
        dmat.dparam$alpha <- setNames(vector(mode = "list", length = nUtype.alpha), Utype.alpha)

        for(iAlpha in Utype.alpha){ ## iAlpha <- Utype.alpha[1]
            dmat.dparam$alpha[[iAlpha]] <- matrix(as.numeric(latent %in% type.alpha[type.alpha$param==iAlpha,"Y"]),
                                                  nrow = n.cluster, ncol = n.latent, byrow = TRUE,
                                                  dimnames = list(NULL,latent))
        }
    }

    if("K" %in% type$detail){
        type.K <- type[type$detail == "K",]
        Utype.K <- unique(type.K$param)
        nUtype.K <- length(Utype.K)

        for(iK in Utype.K){ ## iK <- Utype.K[1]
            iType.K <- type.K[type.K$param==iK,]
            dmat.dparam$K[[iK]] <- matrix(0,
                                          nrow = n.cluster, ncol = n.endogenous,
                                          dimnames = list(NULL,endogenous))

            for(iY in 1:NROW(iType.K)){ ## iY <- 5
                iEndo <- match(iType.K$Y[iY],endogenous)
                dmat.dparam$K[[iK]][,iEndo] <- NA
                dmat.dparam$K[[iK]][X[X$XXendogenousXX==endogenous[iEndo],"XXclusterXX"],iEndo] <- X[X$XXendogenousXX==endogenous[iEndo],iType.K$X[1]]
            }
        }
    }

    if("Gamma" %in% type$detail){
        type.Gamma <- type[type$detail == "Gamma",]
        Utype.Gamma <- unique(type.Gamma$param)
        nUtype.Gamma <- length(Utype.Gamma)
        dmat.dparam$Gamma <- setNames(vector(mode = "list", length = nUtype.Gamma), Utype.Gamma)

        for(iGamma in Utype.Gamma){ ## iGamma <- Utype.Gamma[1]
            iType.Gamma <- type.Gamma[type.Gamma$param==iGamma,]
            dmat.dparam$Gamma[[iGamma]] <- matrix(0,
                                                  nrow = n.cluster, ncol = n.latent,
                                                  dimnames = list(NULL,latent))

            for(iEta in 1:NROW(iType.Gamma)){ ## iLatent <- 5
                iLatent <- match(iType.Gamma$Y[iEta],latent)
                dmat.dparam$Gamma[[iGamma]][,iLatent] <- X[X$XXendogenousXX==latent[iLatent],iType.Gamma$X[1]]
            }
        }
    }

    if("Lambda" %in% type$detail){
        type.Lambda <- type[type$detail == "Lambda",]
        Utype.Lambda <- unique(type.Lambda$param)
        nUtype.Lambda <- length(Utype.Lambda)
        dmat.dparam$Lambda <- setNames(vector(mode = "list", length = nUtype.Lambda), Utype.Lambda)

        for(iLambda in Utype.Lambda){ ## iLambda <- Utype.Lambda[1]
            dmat.dparam$Lambda[[iLambda]] <- matrix(as.double(object$param$Lambda %in% iLambda),
                                                    nrow = n.latent, ncol = n.endogenous,
                                                    dimnames = list(latent,endogenous))
        }
    }

    if("B" %in% type$detail){
        type.B <- type[type$detail == "B",]
        Utype.B <- unique(type.B$param)
        nUtype.B <- length(Utype.B)
        dmat.dparam$B <- setNames(vector(mode = "list", length = nUtype.B), Utype.B)

        for(iB in Utype.B){ ## iB <- Utype.B[1]
            dmat.dparam$B[[iB]] <- matrix(as.double(object$param$B %in% iB),
                                          nrow = n.latent, ncol = n.latent,
                                          dimnames = list(latent,latent))
        }
    }

    if(any(c("Sigma_var","Sigma_cov") %in% type$detail)){
        type.Sigma <- type[type$detail %in% c("Sigma_var","Sigma_cov"),]
        Utype.Sigma <- unique(type.Sigma$param)
        nUtype.Sigma <- length(Utype.Sigma)
        dmat.dparam$Sigma <- setNames(vector(mode = "list", length = nUtype.Sigma), Utype.Sigma)

        for(iSigma in Utype.Sigma){ ## iSigma <- Utype.Sigma[1]
            dmat.dparam$Sigma[[iSigma]] <- matrix(as.double(object$param$Sigma %in% iSigma),
                                                  nrow = n.endogenous, ncol = n.endogenous,
                                                  dimnames = list(endogenous,endogenous))
        }
    }

    if(any(c("sigma2") %in% type$detail)){
        type.sigma2 <- type[type$detail %in% c("sigma2"),]
        Utype.sigma2 <- unique(type.sigma2$param)
        nUtype.sigma2 <- length(Utype.sigma2)
        dmat.dparam$sigma2 <- setNames(vector(mode = "list", length = nUtype.sigma2), Utype.sigma2)

        for(isigma2 in Utype.sigma2){ ## isigma2 <- Utype.sigma2[1]
            dmat.dparam$sigma2[[isigma2]] <- matrix(as.double(object$param$SigmaParam[,,"sigma2"] %in% isigma2),
                                                    nrow = n.endogenous, ncol = n.endogenous,
                                                    dimnames = list(endogenous,endogenous))
        }
    }

    if(any(c("sigma2k") %in% type$detail)){
        type.sigma2k <- type[type$detail %in% c("sigma2k"),]
        Utype.sigma2k <- unique(type.sigma2k$param)
        nUtype.sigma2k <- length(Utype.sigma2k)
        dmat.dparam$sigma2k <- setNames(vector(mode = "list", length = nUtype.sigma2k), Utype.sigma2k)

        for(isigma2k in Utype.sigma2k){ ## isigma2k <- Utype.sigma2k[1]
            dmat.dparam$sigma2k[[isigma2k]] <- array(0, dim = c(n.endogenous,n.endogenous, 2),
                                                     dimnames = list(endogenous,endogenous, c("X","Y")))
            dmat.dparam$sigma2k[[isigma2k]][,,"X"] <- as.double(object$param$SigmaParam[,,"sigma2kX"] %in% isigma2k)
            dmat.dparam$sigma2k[[isigma2k]][,,"Y"] <- as.double(object$param$SigmaParam[,,"sigma2kY"] %in% isigma2k)
        }
    }

    if(any(c("cor") %in% type$detail)){
        type.cor <- type[type$detail %in% c("cor"),]
        Utype.cor <- unique(type.cor$param)
        nUtype.cor <- length(Utype.cor)
        dmat.dparam$cor <- setNames(vector(mode = "list", length = nUtype.cor), Utype.cor)
        
        for(icor in Utype.cor){ ## icor <- Utype.cor[1]
            dmat.dparam$cor[[icor]] <- matrix(as.double(object$param$SigmaParam[,,"cor"] %in% icor),
                                              nrow = n.endogenous, ncol = n.endogenous,
                                              dimnames = list(endogenous,endogenous))
        }
    }

    if(any(c("Psi_var","Psi_cov") %in% type$detail)){
        type.Psi <- type[type$detail %in% c("Psi_var","Psi_cov"),]
        Utype.Psi <- unique(type.Psi$param)
        nUtype.Psi <- length(Utype.Psi)
        dmat.dparam$Psi <- setNames(vector(mode = "list", length = nUtype.Psi), Utype.Psi)

        for(iPsi in Utype.Psi){ ## iB <- Utype.B[1]
            dmat.dparam$Psi[[iPsi]] <- matrix(as.double(object$param$Psi %in% iPsi),
                                              nrow = n.latent, ncol = n.latent,
                                              dimnames = list(latent,latent))
        }
    }

    ## ** Store derivative with respect to the mean/variance
    if(length(dmat.dparam$nu)+length(dmat.dparam$K) > 0){
        dmu.dparam <- c(dmat.dparam$nu, dmat.dparam$K)
    }else{
        dmu.dparam <- list()
    }

    if(length(dmat.dparam$Sigma) > 0){
        dOmega.dparam <- dmat.dparam$Sigma
    }else{
        dOmega.dparam <- list()
    }

    ## ** pairs of parameters to be considered
    grid.param <- list(mean = .combination(object$Uparam.mean, object$Uparam.mean),
                       var = .combination(object$Uparam.var, object$Uparam.var),
                       hybrid = .combination(object$Uparam.mean, object$Uparam.var))


    ## ** export
    return(c(object,
             list(
                 dmu.dparam = dmu.dparam,
                 dOmega.dparam = dOmega.dparam,
                 dmat.dparam = dmat.dparam,
                 grid.dmoment = grid.param
                 ))
           )
}

## * skeletonDtheta2
#' @rdname skeleton
skeletonDtheta2 <- function(object){
    if(lava.options()$debug){cat("skeletonDtheta2 \n")}

    type.param <- object$type[!is.na(object$type$param),,drop=FALSE]
    grid.param <- list()
    
    ## ** identify all combinations of coefficients with second derivative
    grid.param$mean <- list()

    grid.param$mean$alpha.B <- .combinationDF(type.param,
                                        detail1 = "alpha", name1 = "alpha",
                                        detail2 = "B", name2 = "B")

    grid.param$mean$alpha.Lambda <- .combinationDF(type.param,
                                                 detail1 = "alpha", name1 = "alpha",
                                                 detail2 = "Lambda", name2 = "Lambda")

    grid.param$mean$Gamma.B <- .combinationDF(type.param,
                                        detail1 = "Gamma", name1 = "Gamma",
                                        detail2 = "B", name2 = "B")

    grid.param$mean$Gamma.Lambda <- .combinationDF(type.param,
                                             detail1 = "Gamma", name1 = "Gamma",
                                             detail2 = "Lambda", name2 = "Lambda")
    
    grid.param$mean$Lambda.B <- .combinationDF(type.param,
                                        detail1 = "Lambda", name1 = "Lambda",
                                        detail2 = "B", name2 = "B")

    grid.param$mean$B.B <- .combinationDF(type.param,
                                    detail1 = "B", name1 = "B1",
                                    detail2 = "B", name2 = "B2")

    grid.param$var <- list()
    
    grid.param$var$Psi.Lambda <- .combinationDF(type.param,
                                           detail1 = c("Psi_var","Psi_cov"), name1 = "Psi",
                                           detail2 = "Lambda", name2 = "Lambda")

    grid.param$var$Psi.B <- .combinationDF(type.param,
                                      detail1 = c("Psi_var","Psi_cov"), name1 = "Psi",
                                      detail2 = "B", name2 = "B")

    grid.param$var$Lambda.B <- .combinationDF(type.param,
                                             detail1 = "Lambda", name1 = "Lambda",
                                             detail2 = "B", name2 = "B")

    grid.param$var$Lambda.Lambda <- .combinationDF(type.param,
                                              detail1 = "Lambda", name1 = "Lambda1",
                                              detail2 = "Lambda", name2 = "Lambda2")

    grid.param$var$B.B <- .combinationDF(type.param,
                                        detail1 = "B", name1 = "B1",
                                        detail2 = "B", name2 = "B2")

    grid.param$var$sigma2.cor <- .combinationDF(type.param,
                                                detail1 = "sigma2", name1 = "sigma2",
                                                detail2 = "cor", name2 = "cor")

    grid.param$var$sigma2.sigma2k <- .combinationDF(type.param,
                                                    detail1 = "sigma2", name1 = "sigma2",
                                                    detail2 = "sigma2k", name2 = "sigma2k")

    grid.param$var$cor.sigma2k <- .combinationDF(type.param,
                                                 detail1 = "cor", name1 = "cor",
                                                 detail2 = "sigma2k", name2 = "sigma2k")

    grid.param$var$sigma2k.sigma2k <- .combinationDF(type.param,
                                                     detail1 = "sigma2k", name1 = "sigma2k1",
                                                     detail2 = "sigma2k", name2 = "sigma2k2")

    
    grid.param$mean <- grid.param$mean[lengths(grid.param$mean)>0]
    grid.param$var <- grid.param$var[lengths(grid.param$var)>0]
    
    ## ** create d2mu and d2Omega
    if(length(grid.param$n.mean)>0){
        grid.tempo <- lapply(grid.param$mean, function(x){
            if(NROW(x)>0){
                colnames(x) <- c("x","y")
            }
            return(x)
        })
        collapseGrid <- do.call(rbind, grid.tempo)
        name.tempo <- as.character(unique(collapseGrid[[1]]))
        d2mu <- lapply(name.tempo, function(x){
            iIndex <- which(collapseGrid[[1]]==x)
            v <- vector(mode = "list", length(iIndex))
            names(v) <- collapseGrid[[2]][iIndex]
            return(v)
        })
        names(d2mu) <- name.tempo
    }else{
        d2mu <- list()
    }
    
    if(length(grid.param$var)>0){
        grid.tempo <- lapply(grid.param$var, function(x){
            if(NROW(x)>0){
                colnames(x) <- c("x","y")
            }
            return(x)
        })
        collapseGrid <- do.call(rbind, grid.tempo)
        name.tempo <- as.character(unique(collapseGrid[[1]]))
        d2Omega <- lapply(name.tempo, function(x){
            iIndex <- which(collapseGrid[[1]]==x)
            v <- vector(mode = "list", length(iIndex))
            names(v) <- collapseGrid[[2]][iIndex]
            return(v)
        })
        names(d2Omega) <- name.tempo
    }else{
        d2Omega <- list()
    }

    ## ** update grid.dmoment
    object$grid.dmoment$mean$d2.12 <- FALSE
    object$grid.dmoment$mean$d2.21 <- FALSE
    if(length(grid.param$mean)>0){
        for(iType in names(grid.param$mean)){  ## iType <- names(grid.param$mean)[1]
            for(iRow in 1:NROW(grid.param$mean[[iType]])){ ## iRow <- 1
                iName1 <- grid.param$mean[[iType]][iRow,1]
                iName2 <- grid.param$mean[[iType]][iRow,2]
                iIndex1 <- which((object$grid.dmoment$mean$Var1==iName1)*(object$grid.dmoment$mean$Var2==iName2)==1)
                iIndex2 <- which((object$grid.dmoment$mean$Var1==iName2)*(object$grid.dmoment$mean$Var2==iName1)==1)
                if(length(iIndex1)>0){
                    object$grid.dmoment$mean[iIndex1,"d2.12"] <- TRUE
                }else if(length(iIndex2)>0){
                    object$grid.dmoment$mean[iIndex2,"d2.21"] <- TRUE
                }
            }
        }
    }
    
    object$grid.dmoment$var$d2.12 <- FALSE
    object$grid.dmoment$var$d2.21 <- FALSE
    if(length(grid.param$var)>0){
        for(iType in names(grid.param$var)){  ## iType <- names(grid.param$var)[1]
            for(iRow in 1:NROW(grid.param$var[[iType]])){ ## iRow <- 1
                iName1 <- grid.param$var[[iType]][iRow,1]
                iName2 <- grid.param$var[[iType]][iRow,2]
                iIndex1 <- which((object$grid.dmoment$var$Var1==iName1)*(object$grid.dmoment$var$Var2==iName2)==1)
                iIndex2 <- which((object$grid.dmoment$var$Var1==iName2)*(object$grid.dmoment$var$Var2==iName1)==1)
                if(length(iIndex1)>0){
                    object$grid.dmoment$var[iIndex1,"d2.12"] <- TRUE
                }else if(length(iIndex2)>0){
                    object$grid.dmoment$var[iIndex2,"d2.21"] <- TRUE
                }
            }
        }
    }

    ## ** Parameters in dInformation
    ## all but nu and K parameters
    type.Uparam <- type.param[!is.na(type.param$originalLink),]
    if(any(c("Lambda","B") %in% type.Uparam$detail)){
        type.dInformation <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov","sigma2","sigma2k","cor")
    }else{ ## mean-variance model
        type.dInformation <- c("Psi_var","Sigma_var","Psi_cov","Sigma_cov","sigma2","sigma2k","cor")
    }
    name.param.dInformation <- type.Uparam[type.Uparam$detail %in% type.dInformation,"param"]
    
    grid.dInformation <- expand.grid(X = object$Uparam, Y = object$Uparam, Z = name.param.dInformation, stringsAsFactors = FALSE)
    grid.dInformation$level <- apply(grid.dInformation,1,function(iX){paste0(c(sort(iX[1:2]),iX[3]),collapse="")})
    grid.dInformation$duplicated <- duplicated(grid.dInformation$level)

    level2row <- setNames(which(grid.dInformation$duplicated==FALSE),grid.dInformation[grid.dInformation$duplicated==FALSE,"level"])
    grid.dInformation$reference <- level2row[grid.dInformation$level]
    
    
    ## ** Export
    return(c(object,
             list(d2mu.dparam = d2mu,
                  d2Omega.dparam = d2Omega,
                  grid.d2moment = grid.param,
                  name.param.dInformation = name.param.dInformation,
                  grid.dInformation = grid.dInformation))
           )
}



## * helpers
## ** .combinationDF
.combinationDF <- function(data,
                           detail1, detail2,
                           name1, name2){

    detail <- NULL # [:for CRAN check] subset
    
    if(any(detail1 %in% data$detail) && any(detail2 %in% data$detail) ){
        ls.args <- list(subset(data, subset = detail %in% detail1, select = "param", drop = TRUE),
                        subset(data, subset = detail %in% detail2, select = "param", drop = TRUE))
        names(ls.args) <- c(name1,name2)
    
        return(do.call(.combination, args = ls.args))
        
    }else{
        
        return(numeric(0))
        
    }
}


##----------------------------------------------------------------------
### sCorrect-skeleton.R ends here

