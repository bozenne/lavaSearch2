### sCorrect-hessian2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (14:09) 
## Version: 
## Last-Updated: feb 11 2020 (17:23) 
##           By: Brice Ozenne
##     Update #: 58
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - hessian2
#' @title  Extract the Hessian After Small Sample Correction.
#' @description  Extract the hessian from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' @name hessian2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param as.lava [logical] Should the order and the name of the coefficient be the same as those obtained using coef with type = -1.
#' Only relevant for \code{lvmfit} objects.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return An array containing the second derivative of the likelihood relative to each sample (dim 3)
#' and each pair of model coefficients (dim 1,2).
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' eS.lm <- sCorrect(e.lm, ssc = NA, df = "Satterthwaite")
#' hessian2.tempo <- hessian2(eS.lm)
#' colMeans(score.tempo)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' eS.lvm <- sCorrect(e.lvm, df = NA, ssc = NA)
#' score.tempo <- score2(eS.lvm, indiv = TRUE)
#' range(score.tempo-score(e.lvm, indiv = TRUE))
#'
#' @concept small sample inference
#' @export
`hessian2` <-
  function(object, cluster, as.lava) UseMethod("hessian2")

## * hessian2.sCorrect
#' @rdname hessian2
#' @export
hessian2.sCorrect <- function(object, cluster = NULL, as.lava = TRUE){
    
    ## ** define cluster
    if(is.null(cluster)){
        n.cluster <- object$sCorrect$cluster$n.cluster
        cluster.index <- 1:n.cluster
    }else{
        if(!is.numeric(cluster)){
            data <- object$sCorrect$data
            if(length(cluster)==1){                
                if(cluster %in% names(data) == FALSE){
                    stop("Invalid \'cluster\' argument \n",
                         "Could not find variable \"",cluster,"\" in argument \'data\' \n")
                }
                cluster <- data[[cluster]]
            }
            cluster.index <- as.numeric(factor(cluster, levels = unique(cluster)))            
        }else{
            cluster.index <- as.numeric(factor(cluster, levels = unique(cluster)))
        }

        n.cluster <- length(unique(cluster.index))
    }

    ## ** get hessian
    hessian <- object$sCorrect$hessian
    if(is.null(hessian)){return(NULL)}
    if(!is.null(cluster)){ ## aggregate hessian by cluster
        hessianSave <- hessian
        hessian <- array(0, dim = dim(hessian),
                         dimnames = dimnames(hessian))
        for(i in 1:length(cluster)){
            hessian[,,cluster[i]] <- hessian[,,cluster[i]] + hessianSave[,,i]
        }
    }
    
    ## ** export
    if(as.lava == FALSE){
        name.param <- object$sCorrect$name.param
        hessian <- hessian[names(name.param),names(name.param),,drop=FALSE]
        dimnames(hessian) <- list(as.character(name.param),as.character(name.param),NULL)
    }
    if(!is.null(cluster)){
        index2.cluster <- tapply(1:length(cluster),cluster,list)
        attr(hessian,"cluster") <- names(index2.cluster)
    }
    
    return(hessian)
}


## * .hessian2
#' @title Compute the Hessian Matrix From the Conditional Moments
#' @description Compute the Hessian matrix from the conditional moments.
#' @name hessian2-internal
#' 
#' @details \code{calc_hessian} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.hessian2 <- function(dmu, dOmega, d2mu, d2Omega, epsilon, OmegaM1,
                      missing.pattern, unique.pattern, name.pattern,
                      grid.mean, grid.var, grid.hybrid, name.param,
                      leverage, n.cluster){
    if(lava.options()$debug){cat(".hessian2\n")}

    ## ** Prepare
    n.grid.mean <- NROW(grid.mean)
    n.grid.var <- NROW(grid.var)
    n.grid.hybrid <- NROW(grid.hybrid)
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)

    hessian <- array(0, dim = c(n.param, n.param, n.cluster),
                     dimnames = list(name.param,name.param,NULL))
    if(length(dmu)>0){
        index.mean <- 1:n.grid.mean
    }else{
        index.mean <- NULL
    }
    if(length(dOmega)>0){
        index.var <- 1:n.grid.var
    }else{
        index.var <- NULL
    }
    if(length(dmu)>0 && length(dOmega)>0){
        index.hybrid <- 1:n.grid.hybrid
    }else{
        index.hybrid <- NULL
    } 

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        iEpsilon <- epsilon[iIndex,iY,drop=FALSE]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)
        iLeverage <- leverage[iIndex,iY,drop=FALSE]

        ## *** second derivative relative to the mean parameters
        for(iG in index.mean){ # iG <- 2
            iP1 <- grid.mean[iG,1]
            iP2 <- grid.mean[iG,2]
            if(grid.mean[iG,"d2.12"]){
                term1 <- rowSums((id2mu[[iP1]][[iP2]] %*% iOmegaM1) * iEpsilon)
            }else if(grid.mean[iG,"d2.21"]){
                term1 <- rowSums((id2mu[[iP2]][[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1 <- 0
            }
            term2 <- -rowSums((idmu[[iP1]] %*% iOmegaM1) * idmu[[iP2]])
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] + term1 + term2
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
        ## *** second derivative relative to the variance parameters
        for(iG in index.var){ # iG <- 2
            iP1 <- grid.var[iG,1]
            iP2 <- grid.var[iG,2]

            term1a <- - diag(iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1 %*% idOmega[[iP2]])
            term2 <- - rowSums((iEpsilon %*% iOmegaM1 %*% idOmega[[iP2]] %*% iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1) * iEpsilon)
            if(grid.var[iG,"d2.12"]){
                term1b <- diag(iOmegaM1 %*% id2Omega[[iP1]][[iP2]])
                term3 <- 1/2 * rowSums((iEpsilon %*% iOmegaM1 %*% id2Omega[[iP1]][[iP2]] %*% iOmegaM1) * iEpsilon)
            }else if(grid.var[iG,"d2.21"]){
                term1b <- diag(iOmegaM1 %*% id2Omega[[iP2]][[iP1]])
                term3 <- 1/2 * rowSums((iEpsilon %*% iOmegaM1 %*% id2Omega[[iP2]][[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1b <- 0
                term3 <- 0
            }
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] - 1/2 * rowSums( sweep(1-iLeverage, FUN = "*", STATS = term1a + term1b, MARGIN = 2) ) + term2 + term3
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
        
        ## *** second derivative relative to the mean and variance parameters
        for(iG in index.hybrid){ # iG <- 1
            iP1 <- grid.hybrid[iG,1]
            iP2 <- grid.hybrid[iG,2]

            if(!is.null(idmu[[iP1]]) && !is.null(idOmega[[iP2]])){
                term1 <- - rowSums((idmu[[iP1]] %*% iOmegaM1 %*% idOmega[[iP2]] %*% iOmegaM1) * iEpsilon)
            }else{
                term1 <- 0
            }
            if(!is.null(idmu[[iP2]]) && !is.null(idOmega[[iP1]])){
                term2 <- - rowSums((idmu[[iP2]] %*% iOmegaM1 %*% idOmega[[iP1]] %*% iOmegaM1) * iEpsilon)
            }else{
                term2 <- 0
            }
            
            hessian[iP1,iP2,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE] + term1 + term2
            hessian[iP2,iP1,iIndex] <- hessian[iP1,iP2,iIndex,drop=FALSE]
        }
    }

    ## ** export
    return(hessian)
}

## * .subsetList
.subsetList <- function(object, indexRow, indexCol){
    if(length(object)==0){
        return(object)
    }else{    
        return(lapply(object, FUN = function(iL){iL[indexRow,indexCol,drop=FALSE]}))
    }
}
## * .subsetList2
.subsetList2 <- function(object, indexRow, indexCol){
    if(length(object)==0){
        return(object)
    }else{    
        return(lapply(object, FUN = .subsetList, indexRow = indexRow, indexCol = indexCol))
    }
}


######################################################################
### sCorrect-hessian2.R ends here
