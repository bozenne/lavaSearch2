### information2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (14:17) 
## Version: 
## Last-Updated: mar  7 2018 (11:48) 
##           By: Brice Ozenne
##     Update #: 65
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .information2
#' @title Compute the Expected Information Matrix From the Conditional Moments
#' @description Compute the expected information matrix from the conditional moments
#' when there is no missing value.
#'
#' @param n.corrected [numeric >0] the number of independent epsilon.
#' May not match the number of observations when a small sample correction is used.
#'
#' @details \code{.information2} will perform the computation individually when the
#' argument \code{index.Omega} is null.
#' 
#' @keywords internal
.information2 <- function(dmu, dOmega,
                          Omega, OmegaM1, n.corrected,
                          index.Omega, leverage, n.cluster,
                          grid.meanparam, n.grid.meanparam,
                          grid.varparam, n.grid.varparam,
                          name.param, name.meanparam, name.varparam,
                          param2index, n.param){

### ** Prepare
    test.global <- is.null(index.Omega)
    
    Info <- matrix(0, nrow = n.param, ncol = n.param,
                   dimnames = list(name.param,name.param))
    if(length(dmu)>0){
        index.meanparam <- 1:n.grid.meanparam
    }else{
        index.meanparam <- NULL
    }
    if(length(dOmega)>0){
        index.varparam <- 1:n.grid.varparam
    }else{
        index.varparam <- NULL
    } 
    
### ** Global
    if(test.global){
        ## *** Information relative to the mean parameters
        for(iG in index.meanparam){ # iG <- 1
            iP1 <- grid.meanparam[iG,1]
            iP2 <- grid.meanparam[iG,2]

            Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iP1]] %*% OmegaM1 * dmu[[iP2]])
        }

        ## *** Information realtive to the variance parameters
        for(iG in index.varparam){ # iG <- 1
            iP1 <- grid.varparam[iG,1]
            iP2 <- grid.varparam[iG,2]

            iDiag <- diag(OmegaM1 %*% dOmega[[iP1]] %*% OmegaM1 %*% dOmega[[iP2]])
            Info[iP1,iP2] <- Info[iP1,iP2] + 1/2*sum(iDiag*n.corrected)
        }
    }
    
### ** Individual specific
    if(!test.global){
        ## *** Information relative to the mean parameters
        for(iC in 1:n.cluster){ # iC <- 1
            iIndex <- index.Omega[[iC]]

            for(iG in index.varparam){ # iG <- 1
                iP1 <- grid.meanparam[iG,1]
                iP2 <- grid.meanparam[iG,2]

                Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iP1]][iC,iIndex] %*% OmegaM1[iIndex,iIndex,drop=FALSE] * dmu[[iP2]][iC,iIndex])            
            }

            ## *** Information realtive to the variance parameters
            for(iG in index.varparam){ # iG <- 1
                iP1 <- grid.varparam[iG,1]
                iP2 <- grid.varparam[iG,2]

                iDiag <- diag(OmegaM1[iIndex,iIndex,drop=FALSE] %*% dOmega[[iP1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[iIndex,iIndex,drop=FALSE] %*% dOmega[[iP2]][iIndex,iIndex,drop=FALSE])
                Info[iP1,iP2] <- Info[iP1,iP2] + 1/2 * sum(iDiag * (1 - leverage[iC,iIndex]))            
            }
        }        
    }


    ### ** Make Info a symmetric matrix
    Info <- symmetrize(Info, update.upper = NULL)
    
    ### ** export
    return(Info)
}

## * .d2Information
.d2Information <- function(dmu, d2mu, dOmega, d2Omega,
                           Omega, OmegaM1, n.corrected,
                           index.Omega, leverage, n.cluster,
                           name.param, name.3deriv){

    ### ** prepare
    n.param <- length(name.param)
    index.deriv <- match(name.3deriv, name.param)
        
    dInfo <-  array(0,
                    dim = c(n.param, n.param, length(name.3deriv)),
                    dimnames = list(name.param, name.param, name.3deriv))

    ## ** loop
    for(iDeriv in index.deriv){ # iDeriv <- 4
        for(iP1 in 1:n.param){ # iP1 <- 1
            for(iP2 in iP1:n.param){ # iP2 <- 1
                
                iNameD <- name.param[iDeriv]
                iName1 <- name.param[iP1]
                iName2 <- name.param[iP2]
                 ## cat(iNameD," ",iName1,"",iName2,"\n")
                
                ## *** identify relevant terms
                test.Omega1 <- !is.null(dOmega[[iNameD]]) && !is.null(dOmega[[iName1]]) && !is.null(dOmega[[iName2]])
                test.Omega2a <- !is.null(d2Omega[[iNameD]][[iName1]]) && !is.null(dOmega[[iName2]])
                test.Omega2b <- !is.null(d2Omega[[iName1]][[iNameD]]) && !is.null(dOmega[[iName2]])
                test.Omega3a <- !is.null(d2Omega[[iNameD]][[iName2]]) && !is.null(dOmega[[iName1]])
                test.Omega3b <- !is.null(d2Omega[[iName2]][[iNameD]]) && !is.null(dOmega[[iName1]])
                
                test.mu1a <- !is.null(d2mu[[iNameD]][[iName1]]) && !is.null(dmu[[iName2]])
                test.mu1b <- !is.null(d2mu[[iName1]][[iNameD]]) && !is.null(dmu[[iName2]])
                test.mu2a <- !is.null(d2mu[[iNameD]][[iName2]]) && !is.null(dmu[[iName1]])
                test.mu2b <- !is.null(d2mu[[iName2]][[iNameD]]) && !is.null(dmu[[iName1]])
                test.mu3 <- !is.null(dOmega[[iNameD]]) && !is.null(dmu[[iName1]]) && !is.null(dmu[[iName2]])

                if((test.Omega1 + test.Omega2a + test.Omega2b + test.Omega3a + test.Omega3b + test.mu1a + test.mu1b + test.mu2a + test.mu2b + test.mu3) == 0){
                    next
                }
                
                ## *** extract quantities for computations
                if(test.mu1a){
                    d2mu.D1 <- d2mu[[iNameD]][[iName1]]
                }else if(test.mu1b){
                    d2mu.D1 <- d2mu[[iName1]][[iNameD]]
                }
                if(test.mu2a){
                    d2mu.D2 <- d2mu[[iNameD]][[iName2]]
                }else if(test.mu2b){
                    d2mu.D2 <- d2mu[[iName2]][[iNameD]]
                }
                if(!is.null(dOmega[[iNameD]])){
                    OmegaM1.dOmega.D <- OmegaM1 %*% dOmega[[iNameD]]
                }
                if(!is.null(dOmega[[iName1]])){
                    OmegaM1.dOmega.1 <- OmegaM1 %*% dOmega[[iName1]]
                }
                if(!is.null(dOmega[[iName2]])){
                    OmegaM1.dOmega.2 <- OmegaM1 %*% dOmega[[iName2]]
                }
                if(test.Omega2a){
                    d2Omega.D1 <- d2Omega[[iNameD]][[iName1]]
                }else if(test.Omega2b){
                    d2Omega.D1 <- d2Omega[[iName1]][[iNameD]]
                }
                if(test.Omega3a){
                    d2Omega.D2 <- d2Omega[[iNameD]][[iName2]]
                }else{
                    d2Omega.D2 <- d2Omega[[iName2]][[iNameD]]
                }

                ## *** Compute terms
                if(test.Omega1){
                    iDiag1 <- diag(OmegaM1.dOmega.D %*% OmegaM1.dOmega.1 %*% OmegaM1.dOmega.2)
                    iDiag2 <- diag(OmegaM1.dOmega.1 %*% OmegaM1.dOmega.D %*% OmegaM1.dOmega.2)
                    dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * n.corrected + iDiag2 * n.corrected)
                }

                if(test.Omega2a || test.Omega2b){
                    iDiag <- diag(OmegaM1 %*% d2Omega.D1 %*% OmegaM1.dOmega.2)
                    dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * n.corrected)
                }

                if(test.Omega3a || test.Omega3b){
                    iDiag <- diag(OmegaM1.dOmega.1 %*% OmegaM1 %*% d2Omega.D2)
                    dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * n.corrected)
                }

                if(test.mu1a || test.mu1b){
                    dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% OmegaM1 * dmu[[iName2]])
                }

                if(test.mu2a || test.mu2b){
                    dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu[[iName1]] %*% OmegaM1 * d2mu.D2)
                }

                  
                if(test.mu3){
                    dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu[[iName1]] %*% OmegaM1.dOmega.D %*% OmegaM1 * dmu[[iName2]])
                }
            }
        }
        ## *** Symmetrize
        dInfo[,,iNameD] <- symmetrize(dInfo[,,iNameD], update.upper = NULL)
    }

    ### ** export
    return(dInfo)
}


## * .d2InformationIndiv
.d2InformationIndiv <- function(dmu.dtheta, d2mu.dtheta2,
                                dOmega.dtheta, d2Omega.dtheta2,
                                Omega, ls.indexOmega, leverage,
                                n.param, name.param, name.deriv,
                                n.cluster){

### ** prepare
    index.deriv <- match(name.deriv, name.param)
    clusterSpecific <- !is.null(ls.indexOmega)
    iOmega <- chol2inv(chol(Omega))        

    if(!clusterSpecific){
        Omega.tempo <- Omega
        iOmega.tempo <- iOmega

        ## *** small sample correction               
        df.mean <- Reduce("+",leverage)
        iN.cluster <- as.double(n.cluster - diag(df.mean))
    }
    
    ### ** compute the derivative of the information matrix for each coefficients
    dInfo <-  array(0,
                    dim = c(n.param, n.param, length(name.deriv)),
                    dimnames = list(name.param, name.param, name.deriv))
    
    for(iDeriv in index.deriv){ # iDeriv <- 4
        for(iP1 in 1:n.param){ # iP1 <- 1
            for(iP2 in iP1:n.param){ # iP2 <- 1
                
                iNameD <- name.param[iDeriv]
                iName1 <- name.param[iP1]
                iName2 <- name.param[iP2]

                ##cat(iNameD," ",iName1,"",iName2,"\n")
                
                test.Omega1 <- !is.null(dOmega.dtheta[[iNameD]]) && !is.null(dOmega.dtheta[[iName1]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega2a <- !is.null(d2Omega.dtheta2[[iNameD]][[iName1]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega2b <- !is.null(d2Omega.dtheta2[[iName1]][[iNameD]]) && !is.null(dOmega.dtheta[[iName2]])
                test.Omega3a <- !is.null(d2Omega.dtheta2[[iNameD]][[iName2]]) && !is.null(dOmega.dtheta[[iName1]])
                test.Omega3b <- !is.null(d2Omega.dtheta2[[iName2]][[iNameD]]) && !is.null(dOmega.dtheta[[iName1]])
                
                test.mu1a <- !is.null(d2mu.dtheta2[[iNameD]][[iName1]]) && !is.null(dmu.dtheta[[iName2]])
                test.mu1b <- !is.null(d2mu.dtheta2[[iName1]][[iNameD]]) && !is.null(dmu.dtheta[[iName2]])
                test.mu2a <- !is.null(d2mu.dtheta2[[iNameD]][[iName2]]) && !is.null(dmu.dtheta[[iName1]])
                test.mu2b <- !is.null(d2mu.dtheta2[[iName2]][[iNameD]]) && !is.null(dmu.dtheta[[iName1]])
                test.mu3 <- !is.null(dOmega.dtheta[[iNameD]]) && !is.null(dmu.dtheta[[iName1]]) && !is.null(dmu.dtheta[[iName2]])

                if(test.Omega1 + test.Omega2a + test.Omega2b + test.Omega3a + test.Omega3b + test.mu1a + test.mu1b + test.mu2a + test.mu2b + test.mu3 == 0){
                    next
                }
                
                ## *** Individual specific Omega (e.g. presence of missing values)
                if(clusterSpecific){
                    
                    for(iC in 1:n.cluster){

                        ## prepare
                        Omega.tempo <- Omega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        iOmega.tempo <- iOmega[ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        if(!is.null(dmu.dtheta[[iName1]])){
                            dmu.1 <- dmu.dtheta[[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dmu.dtheta[[iName2]])){
                            dmu.2 <- dmu.dtheta[[iName2]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.mu1a){
                            d2mu.D1 <- d2mu.dtheta2[[iNameD]][[iName1]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.mu1b){
                            d2mu.D1 <- d2mu.dtheta2[[iName1]][[iNameD]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.mu1a){
                            d2mu.D <- d2mu.dtheta2[[iNameD]][[iName2]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.mu1b){
                            d2mu.D <- d2mu.dtheta2[[iName2]][[iNameD]][iC,ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iNameD]])){
                            iOmega.dOmega.D <- iOmega.tempo %*% dOmega.dtheta[[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iName1]])){
                            iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(!is.null(dOmega.dtheta[[iName2]])){
                            iOmega.dOmega.2 <- iOmega.tempo %*% dOmega.dtheta[[iName2]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.Omega2a){
                            d2Omega.D1 <- d2Omega.dtheta2[[iNameD]][[iName1]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }else if(test.Omega2b){
                            d2Omega.D1 <- d2Omega.dtheta2[[iName1]][[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }
                        if(test.Omega3a){
                            d2Omega.D2 <- d2Omega.dtheta2[[iNameD]][[iName2]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }else{
                            d2Omega.D2 <- d2Omega.dtheta2[[iName2]][[iNameD]][ls.indexOmega[[iC]],ls.indexOmega[[iC]],drop=FALSE]
                        }

                        ## small sample correction  
                        iW.cluster <- 1 -  diag(leverage[[iC]])
                        
                        ## compute
                        if(test.Omega1){                            
                            iDiag1 <- diag(iOmega.dOmega.D %*% iOmega.dOmega.1 %*% iOmega.dOmega.2)
                            iDiag2 <- diag(iOmega.dOmega.1 %*% iOmega.dOmega.D %*% iOmega.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iW.cluster + iDiag2 * iW.cluster)
                        }
                        
                        if(test.Omega2a || test.Omega2b){
                            iDiag <- diag(iOmega %*% d2Omega.D1 %*% iOmega.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iW.cluster)
                        }

                        if(test.Omega3a || test.Omega3b){
                            iDiag <- diag(iOmega.dOmega.1 %*% iOmega %*% d2Omega.D2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iW.cluster)
                        }

                        if(test.mu1a || test.mu1b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + d2mu.D1 %*% iOmega %*% t(dmu.2)
                        }

                        if(test.mu2a || test.mu2b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + dmu.1 %*% iOmega %*% t(d2mu.D2)
                        }

                        if(test.mu3){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - dmu.1 %*% iOmega.dOmega.D %*% iOmega.tempo %*% t(dmu.2)
                        }
                    
                }
            }
            
            ## *** Same for all individuals
                if(clusterSpecific == FALSE){

                    ## prepare
                    if(!is.null(dmu.dtheta[[iName1]])){
                        dmu.1 <- dmu.dtheta[[iName1]]
                    }
                    if(!is.null(dmu.dtheta[[iName2]])){
                        dmu.2 <- dmu.dtheta[[iName2]]
                    }
                    if(test.mu1a){
                        d2mu.D1 <- d2mu.dtheta2[[iNameD]][[iName1]]
                    }else if(test.mu1b){
                        d2mu.D1 <- d2mu.dtheta2[[iName1]][[iNameD]]
                    }
                    if(test.mu2a){
                        d2mu.D2 <- d2mu.dtheta2[[iNameD]][[iName2]]
                    }else if(test.mu2b){
                        d2mu.D2 <- d2mu.dtheta2[[iName2]][[iNameD]]
                    }
                    if(!is.null(dOmega.dtheta[[iNameD]])){
                        iOmega.dOmega.D <- iOmega.tempo %*% dOmega.dtheta[[iNameD]]
                    }
                    if(!is.null(dOmega.dtheta[[iName1]])){
                        iOmega.dOmega.1 <- iOmega.tempo %*% dOmega.dtheta[[iName1]]
                    }
                    if(!is.null(dOmega.dtheta[[iName2]])){
                        iOmega.dOmega.2 <- iOmega.tempo %*% dOmega.dtheta[[iName2]]
                    }
                    if(test.Omega2a){
                        d2Omega.D1 <- d2Omega.dtheta2[[iNameD]][[iName1]]
                    }else if(test.Omega2b){
                        d2Omega.D1 <- d2Omega.dtheta2[[iName1]][[iNameD]]
                    }
                    if(test.Omega3a){
                        d2Omega.D2 <- d2Omega.dtheta2[[iNameD]][[iName2]]
                    }else{
                        d2Omega.D2 <- d2Omega.dtheta2[[iName2]][[iNameD]]
                    }

                    ## compute
                    if(test.Omega1){
                        iDiag1 <- diag(iOmega.dOmega.D %*% iOmega.dOmega.1 %*% iOmega.dOmega.2)
                        iDiag2 <- diag(iOmega.dOmega.1 %*% iOmega.dOmega.D %*% iOmega.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * iN.cluster + iDiag2 * iN.cluster)
                    }

                    if(test.Omega2a || test.Omega2b){
                        iDiag <- diag(iOmega %*% d2Omega.D1 %*% iOmega.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.cluster)
                    }

                    if(test.Omega3a || test.Omega3b){
                        iDiag <- diag(iOmega.dOmega.1 %*% iOmega %*% d2Omega.D2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * iN.cluster)
                    }

                    if(test.mu1a || test.mu1b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% iOmega * dmu.2)
                    }

                    if(test.mu2a || test.mu2b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu.1 %*% iOmega * d2mu.D2)
                    }

                  
                    if(test.mu3){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu.1 %*% iOmega.dOmega.D %*% iOmega.tempo * dmu.2)
                    }
            }
            
            }
        }
        dInfo[,,iNameD] <- dInfo[,,iNameD] + t(dInfo[,,iNameD]) - diag(diag(dInfo[,,iNameD]))
    }

    ### ** export
    return(dInfo)
}



##----------------------------------------------------------------------
### information2.R ends here
