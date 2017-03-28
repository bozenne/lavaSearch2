### skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (10:35) 
## Version: 
## Last-Updated: mar 28 2017 (17:10) 
##           By: Brice Ozenne
##     Update #: 907
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
#' @param df.param.all [data.frame] output of \code{\link{coefType}} containing the type of each coefficient.
#' @param param2originalLink [named character vector] matching between the name of the coefficient in lava and their label.
#' @param B,alpha.XGamma,Lambda,Psi [matrix] pre-computed matrix.
#' @param OD [list] the pre-computed quantities for the second derivatives. 
#' @param as.lava [logical] should the name of the links be used to name the coefficient?
#' Otherwise uses the labels (when defined) of each coefficient.
#' @param name.endogenous [character vector] name of the endogenous variables
#' @param name.latent [character vector] name of the latent variables
#' @param p [numeric vector, optional] vector of coefficients at which to evaluate the score.
#' @param data [data.frame, optional] data set.
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
`skeleton` <-
    function(object, ...) UseMethod("skeleton")


## * skeleton.lvm
#' @rdname skeleton
skeleton.lvm <- function(object, as.lava,
                         name.endogenous, name.latent,
                         ...){
    detail <- Y <- NULL ## [:for CRAN check] subset

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
    
### ** prepare
    df.param.all  <- coefType(object, as.lava = FALSE)
    if(as.lava){
        param2originalLink <- subset(df.param.all, subset = !is.na(lava), select = c("originalLink", "param"))
        param2originalLink <- stats::setNames(param2originalLink$originalLink, param2originalLink$param)
    }else{
        param2originalLink <- subset(df.param.all, subset = !is.na(lava), select = "param", drop = TRUE)
        param2originalLink <- stats::setNames(param2originalLink, param2originalLink)
    }
    df.param.detail <- subset(df.param.all, subset = !is.na(detail))
    
    skeleton <- list()
    value <- list()
    skeleton$type <- setNames(df.param.detail$detail, df.param.detail$name)
    skeleton$toUpdate <- stats::setNames(c(rep(FALSE,8),TRUE,TRUE,TRUE,TRUE),
                                         c("nu","K","Lambda","Sigma","alpha","Gamma","B","Psi",
                                           "extra","mu","Omega","param"))
   
### ** Measurement model
    
    ## *** nu
    df.param.nu <-  subset(df.param.detail, subset = detail=="nu", select = c("value", "param", "Y", "name"))
    df.param.nu <- df.param.nu[order(df.param.nu$name),]
    value$nu <- stats::setNames(df.param.nu$value,df.param.nu$Y)
    skeleton$nu <- stats::setNames(param2originalLink[df.param.nu$param],df.param.nu$Y)
    skeleton$toUpdate["nu"] <- any(is.na(value$nu))
    
    ## *** X K
    df.param.K <- subset(df.param.detail, subset = detail == "K", select = c("value", "param", "X", "Y"))
    df.param.K <- df.param.K[order(df.param.K$Y),]

    if(NROW(df.param.K)>0){
        value$K <- stats::setNames(lapply(1:n.endogenous, function(iEndogenous){ # iEndogenous <- 1
            subset(df.param.K, subset = Y == name.endogenous[iEndogenous], select = "value", drop = TRUE)
        }), name.endogenous)
            
        skeleton$K <- stats::setNames(lapply(1:n.endogenous, function(iEndogenous){
            param2originalLink[subset(df.param.K, subset = Y == name.endogenous[iEndogenous], select = "param", drop = TRUE)]
        }), name.endogenous)
    
        skeleton$XK <- stats::setNames(lapply(1:n.endogenous, function(iEndogenous){
            subset(df.param.K, subset = Y == name.endogenous[iEndogenous], select = "X", drop = TRUE)
        }), name.endogenous)
        
        skeleton$toUpdate["K"] <- any(unlist(lapply(value$K,is.na)))
    }
    
    ## *** Lambda
    if(n.latent>0){
        ## define matrix
        value$Lambda <- matrix(0,nrow = n.latent, ncol = n.endogenous,
                               dimnames = list(name.latent,name.endogenous))
        skeleton$Lambda <- matrix(as.character(NA),nrow = n.latent, ncol = n.endogenous,
                                  dimnames = list(name.latent,name.endogenous))
        ## update according to the model
        df.param.Lambda <- subset(df.param.detail, subset = detail == "Lambda", select = c("X","Y","param","value","name"))
        df.param.Lambda$index <- match(df.param.Lambda$X, name.latent) + n.latent * (match(df.param.Lambda$Y, name.endogenous)-1)
        df.param.Lambda <- df.param.Lambda[order(df.param.Lambda$name),]

        ## store in the Lambda matrix the name of the coefficient and their pre-computed values
        dfNA.tempo <- subset(df.param.Lambda, subset = is.na(value))
        skeleton$Lambda[dfNA.tempo$index] <- stats::setNames(param2originalLink[dfNA.tempo$param],dfNA.tempo$Y)
        dfNNA.tempo <- subset(df.param.Lambda, subset = !is.na(value))
        value$Lambda[dfNNA.tempo$index] <- stats::setNames(dfNNA.tempo$value,dfNNA.tempo$Y)
        value$Lambda[!is.na(skeleton$Lambda)] <- NA

        skeleton$toUpdate["Lambda"] <- any(is.na(value$Lambda))
    }

    ## *** Sigma    
    ## define matrix
    value$Sigma <- matrix(0,nrow = n.endogenous, ncol = n.endogenous,
                          dimnames = list(name.endogenous,name.endogenous))
    skeleton$Sigma <- matrix(as.character(NA),nrow = n.endogenous, ncol = n.endogenous,
                             dimnames = list(name.endogenous,name.endogenous))
    
    ## update according to the model
    df.param.Sigma <- subset(df.param.detail,
                             subset = detail %in% c("Sigma_var","Sigma_cov"),
                             select = c("X","Y","param","value","name"))
    df.param.Sigma$index <- match(df.param.Sigma$X, name.endogenous) + n.endogenous*(match(df.param.Sigma$Y, name.endogenous)-1)

    dfNA.tempo <- subset(df.param.Sigma, subset = is.na(value))
    skeleton$Sigma[dfNA.tempo$index] <- param2originalLink[dfNA.tempo$param]
    dfNNA.tempo <- subset(df.param.Sigma, subset = !is.na(value))
    value$Sigma[dfNNA.tempo$index] <- dfNNA.tempo$value

    ## symmetrize
    skeleton$Sigma <- symmetrize(skeleton$Sigma, update.upper = TRUE)
    value$Sigma <- symmetrize(value$Sigma, update.upper = TRUE)
    value$Sigma[!is.na(skeleton$Sigma)] <- NA

    skeleton$toUpdate["Sigma"] <- any(is.na(value$Sigma))

### ** Structural model
    if(n.latent>0){
        ## *** alpha 
        df.param.alpha <-  subset(df.param.detail,
                                  subset = detail=="alpha",
                                  select = c("value","param","Y"))
        value$alpha <- stats::setNames(df.param.alpha$value,df.param.alpha$Y)
        skeleton$alpha <- param2originalLink[stats::setNames(df.param.alpha$param,df.param.alpha$Y)]

        skeleton$toUpdate["alpha"] <- any(is.na(value$alpha))
        
        ## *** X Gamma
        df.param.Gamma <- subset(df.param.detail,
                                 subset = detail=="Gamma",
                                 select = c("value","param","X","Y"))
        
        if(NROW(df.param.Gamma)>0){
            
            value$Gamma <- stats::setNames(lapply(1:n.latent, function(iLatent){
                subset(df.param.Gamma, subset = Y==name.latent[iLatent], select = "value", drop = TRUE)
            }), name.latent)
            
            skeleton$Gamma <- stats::setNames(lapply(1:n.latent, function(iLatent){ # iLatent <- 1
                param2originalLink[subset(df.param.Gamma, subset = Y==name.latent[iLatent], select = "param", drop = TRUE)]
            }), name.latent)
    
            skeleton$XGamma <- stats::setNames(lapply(1:n.latent, function(iLatent){
                subset(df.param.Gamma, subset = Y==name.latent[iLatent], select = "X", drop = TRUE)
            }), name.latent)

            skeleton$toUpdate["Gamma"] <- any(unlist(lapply(value$Gamma,is.na)))
        }
        
        ## *** B
        ## define matrix
        value$B <- matrix(0,nrow = n.latent, ncol = n.latent,
                          dimnames = list(name.latent,name.latent))
        skeleton$B <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                             dimnames = list(name.latent,name.latent))

        if(any("B" %in% df.param.all$detail)){
            ## update according to the model
            df.param.B <- subset(df.param.detail,
                                 subset = detail == "B",
                                 select = c("X", "Y", "param", "value", "name"))
            df.param.B$index <- match(df.param.B$X, name.latent) + n.latent*(match(df.param.B$Y, name.latent)-1)
            dfNA.tempo <- subset(df.param.B, subset = is.na(value))            
            skeleton$B[dfNA.tempo$index] <- param2originalLink[dfNA.tempo$param]
            dfNNA.tempo <- subset(df.param.B, subset = is.na(value))
            value$B[dfNNA.tempo$index] <- dfNNA.tempo$value
            value$B[!is.na(skeleton$B)] <- NA

            skeleton$toUpdate["B"] <- any(is.na(value$B))
        }
    
        ## *** Psi    
        ## define matrix
        value$Psi <- matrix(0,nrow = n.latent, ncol = n.latent,
                            dimnames = list(name.latent,name.latent))
        skeleton$Psi <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                               dimnames = list(name.latent,name.latent))

        ## update according to the model
        df.param.Psi <- subset(df.param.all,
                               subset = detail %in% c("Psi_var","Psi_cov"),
                               select = c("X", "Y", "param", "value", "Y", "name"))
                               
        df.param.Psi$index <- match(df.param.Psi$X, name.latent) + n.latent*(match(df.param.Psi$Y, name.latent)-1)

        dfNA.tempo <- subset(df.param.Psi, subset = is.na(value))      
        skeleton$Psi[dfNA.tempo$index] <- param2originalLink[dfNA.tempo$param]
        dfNNA.tempo <- subset(df.param.Psi, subset = is.na(value))
        value$Psi[dfNNA.tempo$index] <- dfNNA.tempo$value

        ## symmetrize
        skeleton$Psi <- symmetrize(skeleton$Psi, update.upper = TRUE)
        value$Psi <- symmetrize(value$Psi, update.upper = TRUE)
        value$Psi[!is.na(skeleton$Psi)] <- NA

        skeleton$toUpdate["Psi"] <- any(is.na(value$Psi))
    }

### ** prepare matrix for updating the variance parameter according to the adjusted Omega
    index.matrix <- data.frame(index = which(upper.tri(skeleton$Sigma, diag = TRUE)),
                               which(upper.tri(skeleton$Sigma, diag = TRUE), arr.ind = TRUE)
                               )
    index.keep <- intersect(which(df.param.all$detail %in% c("Sigma_var","Sigma_cov","Psi_var","Psi_cov")),
                            which(!is.na(df.param.all$lava)))
    name.var <- df.param.all[index.keep,"name"]

    name.rhs <- paste(name.endogenous[index.matrix[,"row"]],
                      lava.options()$symbols[2],
                      name.endogenous[index.matrix[,"col"]],
                      sep = "")
    n.rhs <- length(name.rhs)
    
    A <- matrix(0, nrow = n.rhs, ncol = length(name.var),
                dimnames = list(name.rhs, name.var))
    vec.Sigma <- skeleton$Sigma[index.matrix$index]
    for(i in which(!is.na(vec.Sigma))){
        A[i, vec.Sigma[i]] <- 1
    }

    if(n.latent>0){
        index.Psi <- rbind(cbind(index = which(value$Psi!=0),
                                 which(value$Psi!=0, arr.ind = TRUE)),
                           cbind(index = which(is.na(value$Psi)),
                                 which(is.na(value$Psi), arr.ind = TRUE))
                           )
    }else{
        index.Psi <- NULL        
    }
    adjustMoment <- list(index.matrix = index.matrix,
                            index.Psi = index.Psi,
                            A = A,
                            name.var = name.var,
                            n.rhs = n.rhs)


### ** export
    return(list(skeleton = skeleton,
                value = value,
                df.param = df.param.all,
                adjustMoment = adjustMoment, 
                param2originalLink = param2originalLink)
           )
}


## * skeleton.lvmfit
#' @rdname skeleton
skeleton.lvmfit <- function(object, param, data,
                            name.endogenous, name.latent,
                            ...){
    
    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)
    n.data <- NROW(data)
    
    skeleton <- object$conditionalMoment$skeleton
    toUpdate <- skeleton$toUpdate
    value <- object$conditionalMoment$value
    
### ** Update skeleton with the current values
    ## *** nu
    if(toUpdate["nu"]){
        index.update <- which(!is.na(skeleton$nu))
        value$nu[index.update] <- param[skeleton$nu[index.update]]
    }
    
    ## *** K
    if(toUpdate["K"]){
        for(iY in 1:n.endogenous){ # iY <- 3
            if(length(skeleton$K[[iY]])>0){
                index.update <- which(!is.na(skeleton$K[[iY]]))
                value$K[[iY]][index.update] <- param[skeleton$K[[iY]][index.update]]
            }
        }
    }

    ## *** Lambda
    if(toUpdate["Lambda"]){
        index.update <- which(!is.na(skeleton$Lambda))
        value$Lambda[index.update] <- param[skeleton$Lambda[index.update]]
    }
    
    ## *** Sigma
    if(toUpdate["Sigma"]){
        index.update <- which(!is.na(skeleton$Sigma))
        value$Sigma[index.update] <- param[skeleton$Sigma[index.update]]
    }
    
    ## *** extra
    if(toUpdate["extra"]){
        ## linear predictor (measurement model without latent variable)
        value$nu.XK <- matrix(NA, nrow = n.data, ncol = n.endogenous, byrow = TRUE,
                                       dimnames = list(NULL,name.endogenous))
        for(iY in 1:n.endogenous){ # iY <- 1
            iY2 <- name.endogenous[iY]
            if(length(value$K[[iY2]])>0){
                value$nu.XK[,iY2] <- value$nu[iY2] + data[,skeleton$XK[[iY2]],drop=FALSE] %*% value$K[[iY2]]
            }else{
                value$nu.XK[,iY2] <- value$nu[iY2]
            }
        }
    }
        
### ** Structural model
    if(n.latent>0){
        ## *** alpha
        if(toUpdate["alpha"]){
            index.update <- which(!is.na(skeleton$alpha))
            value$alpha[index.update] <- param[skeleton$alpha[index.update]]
        }
        
        ## *** Gamma
        if(toUpdate["Gamma"]){
            for(iLatent in 1:n.latent){
                if(length(skeleton$Gamma[[iLatent]])>0){
                    index.update <- which(!is.na(skeleton$Gamma[[iLatent]]))
                    value$Gamma[[iLatent]][index.update] <- param[skeleton$Gamma[[iLatent]][index.update]]
                }
            }
        }
        
        ## *** B
        if(toUpdate["B"] && length(skeleton$B)>0){
            index.update <- which(!is.na(skeleton$B))
            value$B[index.update] <- param[skeleton$B[index.update]]
        }
        
        ## *** Psi
        if(toUpdate["Psi"] && length(skeleton$Psi)>0){
            index.update <- which(!is.na(skeleton$Psi))
            value$Psi[index.update] <- param[skeleton$Psi[index.update]]
        }
        
        ## *** extra
        if(toUpdate["extra"]){
            ## linear predictor (latent variable)
            value$alpha.XGamma <- matrix(NA,nrow = n.data, ncol = n.latent, byrow = TRUE,
                                         dimnames = list(NULL,name.latent))
        
            for(iLatent in 1:n.latent){
                iLatent2 <- name.latent[iLatent]
                if(length(value$Gamma[[iLatent2]])>0){
                    value$alpha.XGamma[,iLatent2] <- value$alpha[iLatent2] + data[,skeleton$XGamma[[iLatent2]],drop=FALSE] %*% value$Gamma[[iLatent2]]
                }else{
                    value$alpha.XGamma[,iLatent2] <- value$alpha[iLatent2]
                }
            }
            value$iIB <- solve(diag(1,n.latent,n.latent)-value$B)            
            value$alpha.XGamma.iIB <- value$alpha.XGamma %*% value$iIB

            ## other
            value$iIB.Lambda <-  value$iIB %*% value$Lambda    
            value$Psi.iIB <- value$Psi %*% value$iIB
            value$tLambda.tiIB.Psi.iIB <- t(value$iIB.Lambda) %*% value$Psi.iIB
        }
    }

  
### ** Export
    return(value)
}


## * skeletonDtheta
#' @rdname skeleton
`skeletonDtheta` <-
    function(object, ...) UseMethod("skeletonDtheta")

## * skeletonDtheta.lvm
#' @rdname skeleton
skeletonDtheta.lvm <- function(object, data,
                               df.param.all, param2originalLink,
                               name.endogenous, name.latent, ...){

    factitious <- marginal <- param <- value <- X <- Y <- NULL ## [:for CRAN check] subset

    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)

    df.param <- subset(df.param.all, subset = is.na(value) & marginal == FALSE & factitious == FALSE)
    Utype.by.detail <- tapply(df.param$detail, df.param$param, function(x){length(unique(x))})
    if(any(Utype.by.detail>1)){
        stop("cannot constrain two coefficients of different types to be equal \n")
    }
    name.param <- subset(df.param, subset = !duplicated(param), select = param, drop = TRUE)
    n.param <- length(name.param)

    name.originalLink <- as.character(param2originalLink)

### ** prepare
    n.data <- NROW(data)
    name.data <- colnames(data)
    
    mean.param <- c("nu","K","alpha","Gamma","Lambda","B")
    vcov.param <- c("Sigma_var","Sigma_cov","Psi_var","Psi_cov","Lambda","B")    
    dmu <- list()
    dOmega <- list()
    dLambda <- list()
    dB <- list()
    dPsi <- list()

    toUpdate <- stats::setNames(vector(mode = "logical", n.param),name.originalLink)
    
    ### ** Compute derivative or prepare for the derivative
    for(iName in name.param){ # iName <- name.param[1]

        iName2 <- as.character(param2originalLink[iName])
        iType <- unique(subset(df.param, subset = (param == iName), select = "detail", drop = TRUE))
        iY <- subset(df.param, subset = param %in% iName, select = Y, drop = TRUE)
        iX <- subset(df.param, subset = param %in% iName, select = X, drop = TRUE)

        ## *** derivative regarding the mean        
        if(iType %in% mean.param){            
            if(iType=="nu"){
                dmu[[iName2]] <- matrix(as.numeric(name.endogenous %in% iY),
                                        nrow = n.data, ncol = n.endogenous, byrow = TRUE,
                                        dimnames = list(NULL, name.endogenous))
                toUpdate[iName2] <- FALSE
            }else if(iType=="K"){
                dmu[[iName2]] <- matrix(0, nrow = n.data, ncol = n.endogenous, byrow = TRUE,
                                        dimnames = list(NULL, name.endogenous))
                for(Y.tempo in unique(iY)){                    
                    dmu[[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
                }
                toUpdate[iName2] <- FALSE
            }else if(iType=="alpha"){
                dmu[[iName2]] <- matrix(as.numeric(name.latent %in% unique(iY)), nrow = n.data, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL, name.latent))                
                toUpdate[iName2] <- TRUE
            }else if(iType=="Gamma"){
                dmu[[iName2]] <- matrix(0, nrow = n.data, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL, name.latent))
                for(Y.tempo in unique(iY)){ # Y.tempo <- "eta"
                    dmu[[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
                }
                toUpdate[iName2] <- TRUE
            }
        }
        
        ## *** derivative regarding the residual variance covariance
        if(iType %in% vcov.param){
            
            if(iType=="Sigma_var"){
                dOmega[[iName2]] <- matrix(0,
                                           nrow = n.endogenous, ncol = n.endogenous, byrow = TRUE,
                                           dimnames = list(name.endogenous, name.endogenous))
                dOmega[[iName2]][match(iX, name.endogenous) + (match(iY, name.endogenous) - 1) * n.endogenous] <- 1
                toUpdate[iName2] <- FALSE
            }else if(iType=="Sigma_cov"){
                dOmega[[iName2]] <- matrix(0,
                                           nrow = n.endogenous, ncol = n.endogenous, byrow = TRUE,
                                           dimnames = list(name.endogenous, name.endogenous))
                dOmega[[iName2]][match(iX, name.endogenous) + (match(iY, name.endogenous) - 1) * n.endogenous] <- 1
                dOmega[[iName2]][match(iY, name.endogenous) + (match(iX, name.endogenous) - 1) * n.endogenous] <- 1
                toUpdate[iName2] <- FALSE
            }
            
        }        

        ## *** matrices
        if(iType=="Lambda"){            
            dLambda[[iName2]] <- matrix(0,
                                        nrow = n.latent, ncol = n.endogenous, byrow = TRUE,
                                        dimnames = list(name.latent, name.endogenous))
            dLambda[[iName2]][match(iX, name.latent) + (match(iY, name.endogenous) - 1) * n.latent] <- 1            
            toUpdate[iName2] <- TRUE
        }else if(iType=="B"){
            dB[[iName2]] <- matrix(0,
                                   nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                   dimnames = list(name.latent, name.latent))
            dB[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            toUpdate[iName2] <- TRUE
        }else if(iType=="Psi_var"){
            dPsi[[iName2]] <- matrix(0,
                                     nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                     dimnames = list(name.latent, name.latent))
            dPsi[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            toUpdate[iName2] <- TRUE
        }else if(iType=="Psi_cov"){
            dPsi[[iName2]] <- matrix(0,
                                     nrow = n.latent, ncol = n.latent, byrow = TRUE,
                                     dimnames = list(name.latent, name.latent))
            dPsi[[iName2]][match(iX, name.latent) + (match(iY, name.latent) - 1) * n.latent] <- 1
            dPsi[[iName2]][match(iY, name.latent) + (match(iX, name.latent) - 1) * n.latent] <- 1            
            toUpdate[iName2] <- TRUE
        } 
    }
### ** export
    return(list(
        dmu = dmu,
        dOmega = dOmega,
        dLambda = dLambda,
        dB = dB,
        dPsi = dPsi,
        toUpdate = toUpdate
    ))
}


## * skeletonDtheta.lvmfit
#' @rdname skeleton
skeletonDtheta.lvmfit <- function(object, name.endogenous, name.latent, ...){

### ** Import information
    n.endogenous <- length(name.endogenous)
    n.latent <- length(name.latent)

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
    
### ** Update partial derivatives

    ## *** mean coefficients
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

    ## *** variance-covariance coefficients
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
        
    

### ** Export
    return(list(dmu = dmu, dOmega = dOmega))

}

## * skeletonDtheta2
#' @rdname skeleton
`skeletonDtheta2` <-
    function(object, ...) UseMethod("skeletonDtheta2")

## * skeletonDtheta2.lvm
#' @rdname skeleton
skeletonDtheta2.lvm <- function(object, data, df.param.all,
                                param2originalLink, name.latent, ...){

    detail <- factitious <- marginal <- param <- value <- Y <- NULL ## [:for CRAN check] subset
    
    df.param <- subset(df.param.all, is.na(value) & marginal == FALSE & factitious == FALSE)
    dfred.param <- subset(df.param, subset = !duplicated(param))
    
    n.latent <- length(name.latent)
    n.data <- NROW(data)

### ** identify all combinations of coefficients with second derivative
    grid.mean <- list()

    grid.mean$alpha.B <- .combinationDF(dfred.param,
                                        detail1 = "alpha", name1 = "alpha",
                                        detail2 = "B", name2 = "B")

    grid.mean$alpha.Lambda <- .combinationDF(dfred.param,
                                             detail1 = "alpha", name1 = "alpha",
                                             detail2 = "Lambda", name2 = "Lambda")

    grid.mean$Gamma.B <- .combinationDF(dfred.param,
                                        detail1 = "Gamma", name1 = "Gamma",
                                        detail2 = "B", name2 = "B")

    grid.mean$Gamma.Lambda <- .combinationDF(dfred.param,
                                             detail1 = "Gamma", name1 = "Gamma",
                                             detail2 = "Lambda", name2 = "Lambda")
    
    grid.mean$Lambda.B <- .combinationDF(dfred.param,
                                        detail1 = "Lambda", name1 = "Lambda",
                                        detail2 = "B", name2 = "B")

    grid.mean$B.B <- .combinationDF(dfred.param,
                                    detail1 = "B", name1 = "B1",
                                    detail2 = "B", name2 = "B2")

    n.mean <- lapply(grid.mean, NROW)
    

    grid.vcov <- list()
    
    grid.vcov$Psi.Lambda <- .combinationDF(dfred.param,
                                           detail1 = c("Psi_var","Psi_cov"), name1 = "Psi",
                                           detail2 = "Lambda", name2 = "Lambda")

    grid.vcov$Psi.B <- .combinationDF(dfred.param,
                                      detail1 = c("Psi_var","Psi_cov"), name1 = "Psi",
                                      detail2 = "B", name2 = "B")

    grid.vcov$Lambda.B <- .combinationDF(dfred.param,
                                         detail1 = "Lambda", name1 = "Lambda",
                                         detail2 = "B", name2 = "B")

    grid.vcov$Lambda.Lambda <- .combinationDF(dfred.param,
                                              detail1 = "Lambda", name1 = "Lambda1",
                                              detail2 = "Lambda", name2 = "Lambda2")

    grid.vcov$B.B <- .combinationDF(dfred.param,
                                    detail1 = "B", name1 = "B1",
                                    detail2 = "B", name2 = "B2")
    
    n.vcov <- lapply(grid.vcov, NROW)
    
### ** convert back to lava names
    grid.mean <- lapply(grid.mean, function(x){ ## x <- grid.mean[[2]]
        if(length(x)>0){
            x[,1] <- param2originalLink[x[,1]]
            x[,2] <- param2originalLink[x[,2]]
        }
        return(x)
    })

    grid.vcov <- lapply(grid.vcov, function(x){ ## x <- grid.vcov[[2]]
        if(length(x)>0){
            x[,1] <- param2originalLink[x[,1]]
            x[,2] <- param2originalLink[x[,2]]
        }
        return(x)
    })

### ** prepare export
    if(any(unlist(n.mean)>0)){
        xx <- lapply(grid.mean, function(x){
            if(NROW(x)>0){
                colnames(x) <- c("x","y")
            }
            return(x)
        })
        collapseGrid <- do.call(rbind, xx)
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
    
    if(any(unlist(n.vcov)>0)){
        xx <- lapply(grid.vcov, function(x){
            if(NROW(x)>0){
                colnames(x) <- c("x","y")
            }
            return(x)
        })
        collapseGrid <- do.call(rbind, xx)
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
    
    ## ** prepare alpha.B and alpha.Lambda
    if(any(df.param$detail == "alpha")){
        name.alpha <- subset(df.param, subset = !duplicated(param) & detail == "alpha", select = "param", drop = TRUE)
        ls.Malpha <- list()
        for(iName in name.alpha){ # iName <- name.alpha[1]

            iParam <- df.param[df.param$name == iName, "param"]
            iY <- subset(df.param, subset = param %in% iParam, select = Y, drop = TRUE)
            ls.Malpha[[iName]] <- matrix(as.numeric(name.latent %in% unique(iY)),
                                         nrow = n.data, ncol = n.latent, byrow = TRUE,
                                         dimnames = list(NULL, name.latent))
            
        }
    }
    if(n.mean$alpha.B>0){
        for(iP in 1:n.mean$alpha.B){ ## iP <- 1
            iName1 <- grid.mean$alpha.B[iP,"alpha"]
            iName2 <- grid.mean$alpha.B[iP,"B"]
            
            d2mu[[iName1]][[iName2]] <- ls.Malpha[[iName1]]
        }
    }
    if(n.mean$alpha.Lambda>0){
        for(iP in 1:n.mean$alpha.Lambda){ ## iP <- 1
            iName1 <- grid.mean$alpha.Lambda[iP,"alpha"]
            iName2 <- grid.mean$alpha.Lambda[iP,"Lambda"]
            
            d2mu[[iName1]][[iName2]] <- ls.Malpha[[iName1]]
        }
    }
    
    ## ** Store X for Gamma
    if(n.mean$Gamma.Lambda>0){
        for(iP in 1:n.mean$Gamma.Lambda){ ## iP <- 1
            iName1 <- grid.mean$Gamma.Lambda[iP,"Gamma"]
            iName2 <- grid.mean$Gamma.Lambda[iP,"Lambda"]

            iParam <- df.param[df.param$name == iName1, "param"]
            iX <- subset(df.param.all, subset = param %in% iParam, select = "X", drop = TRUE)
            iY <- subset(df.param.all, subset = param %in% iParam, select = "Y", drop = TRUE)

            d2mu[[iName1]][[iName2]] <- matrix(0, nrow = n.data, ncol = n.latent, byrow = TRUE,
                                               dimnames = list(NULL, name.latent))
            for(Y.tempo in unique(iY)){
                d2mu[[iName1]][[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
            }
        }
    }
    
    if(n.mean$Gamma.B>0){
        for(iP in 1:n.mean$Gamma.B){ ## iP <- 1
            iName1 <- grid.mean$Gamma.B[iP,"Gamma"]
            iName2 <- grid.mean$Gamma.B[iP,"B"]
            
            iParam <- df.param[df.param$name == iName1, "param"]
            iX <- subset(df.param.all, subset = param %in% iParam, select = "X", drop = TRUE)
            iY <- subset(df.param.all, subset = param %in% iParam, select = "Y", drop = TRUE)

            d2mu[[iName1]][[iName2]] <- matrix(0, nrow = n.data, ncol = n.latent, byrow = TRUE,
                                               dimnames = list(NULL, name.latent))
            for(Y.tempo in unique(iY)){
                d2mu[[iName1]][[iName2]][,Y.tempo] <- rowSums(data[,iX[iY == Y.tempo],drop=FALSE])
            }
        }
    }

### ** export
    toUpdate <- unlist(lapply(c(n.mean,n.vcov), function(x){x>0}))
    return(list(grid.mean = grid.mean,
                n.mean = n.mean,                
                grid.vcov = grid.vcov,
                n.vcov = n.vcov,
                d2mu = d2mu,
                d2Omega = d2Omega,
                toUpdate = toUpdate
                ))
}

## * skeletonDtheta2.lvmfit
#' @rdname skeleton
skeletonDtheta2.lvmfit <- function(object, name.endogenous, name.latent, ...){
    
### ** Import information
    n.endogenous <- length(name.endogenous)

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
    
### ** second order partial derivatives
    if(any(toUpdate)){
        
        ## *** mean coefficients        
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

        ## *** variance-covariance coefficients
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

### ** Export
    return(list(d2mu = d2mu, d2Omega = d2Omega))

}

## * .combination
#' @title Form all Unique Combinations Between two Vectors
#' @description Form all unique combinations between two vectors (removing symmetric combinations).
#' @name combination
#'
#' @param ... [vectors] elements to be combined.
#'
#' @return A matrix, each row being a different combination.
#' 
#' @examples
#' .combination <- lavaSearch2:::.combination
#' 
#' .combination(1,1)
#' .combination(1:2,1:2)
#' .combination(c(1:2,1:2),1:2)
#' 
#' .combination(alpha = 1:2, beta = 3:4)
#'
#' @keywords internal
.combination <- function(...){

    ## ** normalize arguments
    dots <- list(...)
    if(length(dots)!=2){
        stop("can only handle two vectors \n")
    }
    test.null <- unlist(lapply(dots,is.null))    
    if(any(test.null)){
        return(NULL)
    }
    dots <- lapply(dots,unique)

    ## ** form all combinations
    grid <- expand.grid(dots, stringsAsFactors = FALSE) 
    
    ## ** remove combinations (b,a) when (a,b) is already there
    name1 <- paste0(grid[,1],grid[,2])
    name2 <- paste0(grid[,2],grid[,1])

    if(NROW(grid)>1 && any(name1 %in% name2)){ 

        n.grid <- NROW(grid)
        test.duplicated <- c(FALSE,sapply(2:n.grid, function(iG){
            any(name2[iG] %in% name1[1:(iG-1)]) ## find duplicates
        }))

        grid <- grid[test.duplicated==FALSE,]
    }

    ## ** export
    return(grid)        
}


## * .combinationDF
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
### skeleton.R ends here

