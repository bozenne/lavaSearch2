### skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (10:35) 
## Version: 
## Last-Updated: nov  8 2017 (14:57) 
##           By: Brice Ozenne
##     Update #: 112
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


## * skeleton.lvmfit
#' @rdname skeleton
#' @export
skeleton.lvmfit <- function(object, as.lava, update.value,
                            p, data,
                            name.endogenous, n.endogenous,
                            name.latent, n.latent){

### ** prepare
    dt.param.all <- coefType(object, as.lava = FALSE)

    if(as.lava){
        dt.convert <- dt.param.all[is.na(value)][!duplicated(param),.(param,lava,originalLink)]
        param2originalLink <- dt.convert[!is.na(lava)][,setNames(originalLink,param)]
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

### ** update value with the current values
    if(update.value){

        value <- .updateValueWithSkeleton(param = p, data = data,
                                          n.latent = n.latent, name.latent = name.latent,
                                          n = NROW(data), n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                          skeleton = skeleton, value = value)
        
    }
    
### ** export
    return(list(skeleton = skeleton,
                value = value,
                dt.param = dt.param.all,
                param2originalLink = param2originalLink)
           )
}

## * .updateValueWithSkeleton
.updateValueWithSkeleton <- function(param, data,
                                     n.latent, name.latent,
                                     n, n.endogenous, name.endogenous,
                                     skeleton, value){

### ** Measurement model
    ## *** nu
    index.update <- which(!is.na(skeleton$nu))
    value$nu[index.update] <- param[skeleton$nu[index.update]]

    ## *** K 
    for(iY in 1:n.endogenous){ # iY <- 3
        if(length(skeleton$K[[iY]])>0){
            index.update <- which(!is.na(skeleton$K[[iY]]))
            value$K[[iY]][index.update] <- param[skeleton$K[[iY]][index.update]]
        }
    }

    ## *** Lambda
    if(n.latent>0){
        index.update <- which(!is.na(skeleton$Lambda))
        value$Lambda[index.update] <- param[skeleton$Lambda[index.update]]
    }
    
    ## *** Sigma
    index.update <- which(!is.na(skeleton$Sigma))
    value$Sigma[index.update] <- param[skeleton$Sigma[index.update]]

    ## *** linear predictor
    value$nu.XK <- matrix(NA,nrow = n, ncol = n.endogenous, byrow = TRUE,
                          dimnames = list(NULL,name.endogenous))
    for(iY in 1:n.endogenous){ # iY <- 1
        if(length(value$K[[iY]])>0){
            value$nu.XK[,iY] <- value$nu[iY] + data[,skeleton$XK[[iY]],drop=FALSE] %*% value$K[[iY]]
        }else{
            value$nu.XK[,iY] <- value$nu[iY]
        }
    }
        
### ** Structural model
    if(n.latent>0){
        ## *** alpha
        index.update <- which(!is.na(skeleton$alpha))
        value$alpha[index.update] <- param[skeleton$alpha[index.update]]

        ## *** Gamma
        for(iLatent in 1:n.latent){
            if(length(skeleton$Gamma[[iLatent]])>0){
                index.update <- which(!is.na(skeleton$Gamma[[iLatent]]))
                value$Gamma[[iLatent]][index.update] <- param[skeleton$Gamma[[iLatent]][index.update]]
            }
        }
        
        ## *** B
        index.update <- which(!is.na(skeleton$B))
        value$B[index.update] <- param[skeleton$B[index.update]]

        ## *** Psi
        index.update <- which(!is.na(skeleton$Psi))
        value$Psi[index.update] <- param[skeleton$Psi[index.update]]

        ## *** linear predictor
        value$alpha.XGamma <- matrix(NA,nrow = n, ncol = n.latent, byrow = TRUE,
                                     dimnames = list(NULL,name.latent))
        for(iLatent in 1:n.latent){

            if(length(value$Gamma[[iLatent]])>0){
                value$alpha.XGamma[,iLatent] <- value$alpha[iLatent] + data[,skeleton$XGamma[[iLatent]],drop=FALSE] %*% value$Gamma[[iLatent]]
            }else{
                value$alpha.XGamma[,iLatent] <- value$alpha[iLatent]
            }
        }
    }

### ** Export
    return(value)

}



##----------------------------------------------------------------------
### skeleton.R ends here
