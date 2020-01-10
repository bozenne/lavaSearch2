### sCorrect-information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (14:17) 
## Version: 
## Last-Updated: jan 10 2020 (13:24) 
##           By: Brice Ozenne
##     Update #: 364
##----------------------------------------------------------------------
## 
### Commentary: 
## Compute information, hessian, and first derivative of information
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - information2
#' @title  Extract The Full Information Matrix After Small Sample Correction. 
#' @description  Extract the full information matrix from  from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{lava::information} but with small sample correction.
#' @name information2
#'
#' @param object a linear model or a latent variable model
#' @param ... arguments to be passed to \code{vcov2}.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix.
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
#' info.tempo <- vcov2(e.lm, bias.correct = TRUE)
#' info.tempo[names(coef(e.lm)),names(coef(e.lm))] - vcov(e.lm)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov.tempo <- vcov2(e.lvm, bias.correct = FALSE)
#' round(vcov.tempo %*% information(e.lvm), 5)
#'
#' @concept small sample inference
#' @export
`information2` <-
  function(object, param, data, ssc) UseMethod("information2")

## * information2.lm
#' @rdname information2
#' @export
information2.lm <- function(object, param = NULL, data = NULL,
                            ssc = lava.options()$ssc){

    if(is.null(object$sCorrect) || !is.null(param) || !is.null(data) || !identical(object$sCorrect$ssc$type,ssc)){
        object <- sCorrect(object, param = param, data = data, ssc = ssc, df = NA)
    }

    return(object$sCorrect$information)
}

## * information2.gls
#' @rdname information2
#' @export
information2.gls <- information2.lm

## * information2.lme
#' @rdname information2
#' @export
information2.lme <- information2.lm

## * information2.lvmfit
#' @rdname information2
#' @export
information2.lvmfit <- information2.lm

## * information2.sCorrect
#' @rdname information2
information2.sCorrect <- function(object, param = NULL, data = NULL,
                                 ssc = object$sCorrect$ssc$type){
    class(object) <- setdiff(class(object),"sCorrect")
    return(information2(object, param = param, data = data, ssc = ssc))
}

## * information.sCorrect
#' @rdname information2
information.sCorrect <- information2.sCorrect


## * .information2
#' @title Compute the Expected Information Matrix From the Conditional Moments
#' @description Compute the expected information matrix from the conditional moments.
#' @name information2-internal
#' 
#' @details \code{calc_information} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.information2 <- function(dmu, dOmega, OmegaM1,
                          missing.pattern, unique.pattern, name.pattern,
                          grid.mean, grid.var, name.param,
                          leverage, n.cluster){
    if(lava.options()$debug){cat(".information2\n")}
    
    ## ** Prepare
    n.grid.mean <- NROW(grid.mean)
    n.grid.var <- NROW(grid.var)
    n.param <- length(name.param)
    n.pattern <- length(name.pattern)

    Info <- matrix(0, nrow = n.param, ncol = n.param,
                   dimnames = list(name.param,name.param))
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

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iN.corrected <- n.cluster - colSums(leverage[iIndex,iY,drop=FALSE])

        ## *** Information relative to the mean parameters
        for(iG in index.mean){ # iG <- 1
            iP1 <- grid.mean[iG,1]
            iP2 <- grid.mean[iG,2]

            Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iP1]][iIndex,iY,drop=FALSE] %*% iOmegaM1 * dmu[[iP2]][iIndex,iY,drop=FALSE])
        }

        ## *** Information realtive to the variance parameters
        for(iG in index.var){ # iG <- 1
            iP1 <- grid.var[iG,1]
            iP2 <- grid.var[iG,2]
            iDiag <- diag(iOmegaM1 %*% dOmega[[iP1]][iY,iY,drop=FALSE] %*% iOmegaM1 %*% dOmega[[iP2]][iY,iY,drop=FALSE])
            Info[iP1,iP2] <- Info[iP1,iP2] + 1/2*sum(iDiag*iN.corrected)
        }        
    }

    ## ** Make Info a symmetric matrix
    Info <- symmetrize(Info, update.upper = NULL)
        
    ## ** export
    return(Info)
}





