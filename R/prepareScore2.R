### prepareScore2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: nov  6 2017 (11:23) 
##           By: Brice Ozenne
##     Update #: 160
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @name prepareScore2
#' @export
`prepareScore2` <-
  function(x, ...) UseMethod("prepareScore2")

#' @rdname prepareScore2
#' @export
prepareScore2.lvmfit <- function(x, data = NULL){

### ** extract parameters
    dt.param.all <- coefType(x, as.lava = FALSE)
    dt.param <- dt.param.all[is.na(value) & factice == FALSE]
    Utype.by.detail <- dt.param[,.(n.type = length(unique(detail))), by = param][["n.type"]]
    if(any(Utype.by.detail>1)){
        stop("cannot constrain two parameters of different types to be equal \n")
    }

### ** extract data
    if(is.null(data)){
        data <- model.frame(x)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    
### ** prepare
    name.param <- dt.param[!duplicated(param),param]
    param2originalLink <- dt.param.all[!is.na(param)&!is.na(lava)][,setNames(originalLink,param)]
    name.originalLink <- as.character(param2originalLink)
    n.param <- length(name.param)
    
    n <- NROW(data)

    name.endogenous <- endogenous(x)
    n.endogenous <- length(name.endogenous)
    name.latent <- latent(x)
    n.latent <- length(name.latent)
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
  
### ** prepare the computation of the residual variance covariance matrix
    skeleton <- list()
    
    ## Sigma
    Sigma <- matrix(0,nrow = n.endogenous, ncol = n.endogenous,
                    dimnames = list(name.endogenous,name.endogenous))
    skeleton$Sigma <- matrix(as.character(NA),nrow = n.endogenous, ncol = n.endogenous,
                             dimnames = list(name.endogenous,name.endogenous))
    dt.tempo <- dt.param.all[detail %in% c("Sigma_var","Sigma_cov"), .(index = which(name.endogenous %in% X) + n.endogenous*(which(name.endogenous %in% Y)-1), param, value), by = name]
    skeleton$Sigma[dt.tempo[is.na(value)]$index] <- param2originalLink[dt.tempo[is.na(value),param]]
    Sigma[dt.tempo[!is.na(value)]$index] <- dt.tempo[!is.na(value),value]

    Sigma[upper.tri(Sigma, diag = FALSE)] <- Sigma[lower.tri(Sigma, diag = FALSE)]
    skeleton$Sigma[upper.tri(Sigma, diag = FALSE)] <- skeleton$Sigma[lower.tri(Sigma, diag = FALSE)]
    Sigma[!is.na(skeleton$Sigma)] <- NA
    
    ## Psi
    Psi <- matrix(0,nrow = n.latent, ncol = n.latent,
                  dimnames = list(name.latent,name.latent))
    skeleton$Psi <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                           dimnames = list(name.latent,name.latent))
    if(n.latent>0){
        dt.tempo <- dt.param.all[detail %in% c("Psi_var","Psi_cov"), .(index = which(name.latent %in% X) + n.latent*(which(name.latent %in% Y)-1), param, value), by = name]
        skeleton$Psi[dt.tempo[is.na(value)]$index] <- param2originalLink[dt.tempo[is.na(value),param]]
        Psi[dt.tempo[!is.na(value)]$index] <- dt.tempo[!is.na(value),value]

        Psi[upper.tri(Psi, diag = FALSE)] <- Psi[lower.tri(Psi, diag = FALSE)]
        skeleton$Psi[upper.tri(Psi, diag = FALSE)] <- skeleton$Psi[lower.tri(Psi, diag = FALSE)]
        Psi[!is.na(skeleton$Psi)] <- NA
    }

    ## B
    B <- matrix(0,nrow = n.latent, ncol = n.latent,
                dimnames = list(name.latent,name.latent))
    skeleton$B <- matrix(as.character(NA),nrow = n.latent, ncol = n.latent,
                         dimnames = list(name.latent,name.latent))
    if(n.latent>0){
        if(any("B" %in% dt.param.all$detail)){
            dt.tempo <- dt.param.all[detail == "B", .(index = which(name.latent %in% X) + n.latent*(which(name.latent %in% Y)-1), param, value), by = name]
            skeleton$B[dt.tempo[is.na(value)]$index] <- param2originalLink[dt.tempo[is.na(value),param]]
            B[dt.tempo[!is.na(value)]$index] <- dt.tempo[!is.na(value),value]
            B[!is.na(skeleton$B)] <- NA            
        }
    }

    ## Lambda
    Lambda <- matrix(0,nrow = n.latent, ncol = n.endogenous,
                     dimnames = list(name.latent,name.endogenous))
    skeleton$Lambda <- matrix(as.character(NA),nrow = n.latent, ncol = n.endogenous,
                              dimnames = list(name.latent,name.endogenous))
    if(n.latent>0){
        dt.tempo <- dt.param.all[detail %in% "Lambda", .(index = which(name.latent %in% X) + n.latent*(which(name.endogenous %in% Y)-1), param, value), by = name]
        skeleton$Lambda[dt.tempo[is.na(value)]$index] <- param2originalLink[dt.tempo[is.na(value),param]]
        Lambda[dt.tempo[!is.na(value)]$index] <- dt.tempo[!is.na(value),value]
        Lambda[!is.na(skeleton$Lambda)] <- NA
    }

    ## alpha + X Gamma
    alpha <- dt.param.all[detail=="alpha",value,by="name"][["value"]]
    skeleton$alpha <- param2originalLink[dt.param.all[detail=="alpha",value,by="param"][["param"]]]

    dt.tempo <- dt.param.all[detail=="Gamma",value,by="Y"]
    Gamma <- setNames(lapply(1:n.latent, function(iLatent){
        dt.tempo$value[dt.tempo$Y==name.latent[iLatent]]
    }), name.latent)


    dt.tempo <- dt.param.all[detail=="Gamma",param,by="Y"]
    skeleton$Gamma <- setNames(lapply(1:n.latent, function(iLatent){
        param2originalLink[dt.tempo$param[dt.tempo$Y==name.latent[iLatent]]]
    }), name.latent)
    
    dt.tempo <- dt.param.all[detail=="Gamma",X,by="Y"]
    skeleton$XGamma <- setNames(lapply(1:n.latent, function(iLatent){
        dt.tempo$X[dt.tempo$Y==name.latent[iLatent]]
    }), name.latent) 

### ** compute derivative or prepare for the derivative
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
        skeleton = skeleton,
        alpha = alpha,
        Gamma = Gamma,
        Sigma = Sigma,
        Psi = Psi,
        B = B,
        Lambda = Lambda,
        type = type,
        toUpdate = toUpdate
    ))
}

#----------------------------------------------------------------------
### prepareScore2.R ends here
