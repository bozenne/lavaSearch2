### sCorrect-getVarCov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (11:00) 
## Version: 
## Last-Updated: nov 18 2019 (15:37) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getVarCov2
#' @title Reconstruct the Conditional Variance Covariance Matrix After Small Sample Correction.
#' @description Reconstruct the conditional variance covariance matrix from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects. 
#' For \code{gls} and \code{lme} object, it is only compatible with specific correlation and variance structure.
#' It is similar to \code{nlme::getVarCov} but with small sample correction (if any).
#' 
#' @name getVarCov2
#'
#' @param object a \code{gls} or \code{lme} object
#' @param param [numeric vector] values for the model parameters.
#' @param data [data.frame] the data set.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param ... [internal] only used by the generic method.
#' 
#' @details The compound symmetry variance-covariance matrix in a gls model is of the form:
#' \tabular{cccc}{
#' \eqn{\Sigma =} \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \rho} \tab \eqn{\sigma^2 \rho} \cr
#' \tab . \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \rho} \cr
#' \tab . \tab . \tab \eqn{\sigma^2}
#' }
#'
#' The unstructured variance-covariance matrix in a gls model is of the form:
#'  \tabular{cccc}{
#' \eqn{\Sigma =} \tab \eqn{\sigma^2} \tab \eqn{\sigma^2 \sigma_2 \rho_{1,2}} \tab \eqn{\sigma^2 \sigma_3 \rho_{1,3}} \cr
#' \tab . \tab \eqn{\sigma^2 \sigma_2^2} \tab \eqn{\sigma^2 \sigma_2 \sigma_3 \rho_{2,3}} \cr
#' \tab . \tab . \tab \eqn{\sigma^2 \sigma_3^2}
#' }
#' @return A list containing the residual variance-covariance matrix in the element Omega.
#' 
#' @examples
#' 
#' ## simulate data 
#' library(nlme)
#' n <- 5e1
#' mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
#' latent(mSim) <- ~eta
#' transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
#' set.seed(10)
#' dW <- lava::sim(mSim,n,latent = FALSE)
#' dW <- dW[order(dW$Id),,drop=FALSE]
#' dL <- reshape2::melt(dW,id.vars = c("G","Id"), variable.name = "time")
#' dL <- dL[order(dL$Id),,drop=FALSE]
#' dL$Z1 <- rnorm(NROW(dL))
#' dL$time.num <- as.numeric(as.factor(dL$time))
#'
#' #### linear regression ####
#' e1.lm <- lm(Y1 ~ G, data = dW)
#' getVarCov2(e1.lm)
#' 
#' #### iid model #### 
#' e1.gls <- nlme::gls(Y1 ~ G, data = dW, method = "ML")
#' getVarCov2(e1.gls, cluster = 1:n)$Omega
#' 
#' #### heteroschedasticity ####
#' dW$group <- rbinom(n, size = 1, prob = 1/2)
#' dW$repetition <- as.numeric(as.factor(dW$group))
#' e2a.gls <- nlme::gls(Y1 ~ G, data = dW, method = "ML",
#'                     weights = varIdent(form =~ repetition|group))
#' getVarCov2(e2a.gls, cluster = 1:n)$Omega
#'
#' 
#' e2b.gls <- nlme::gls(value ~ 0+time + time:G,
#'                    weight = varIdent(form = ~ time.num|time),
#'                    data = dL, method = "ML")
#' getVarCov2(e2b.gls, cluster = "Id")$Omega
#'
#' #### compound symmetry ####
#' e3.gls <- nlme::gls(value ~ time + G,
#'                    correlation = corCompSymm(form = ~1| Id),
#'                    data = dL, method = "ML")
#' getVarCov2(e3.gls)$Omega
#' 
#' #### unstructured ####
#' e4.gls <- nlme::gls(value ~ time,
#'                     correlation = corSymm(form = ~time.num| Id),
#'                     weight = varIdent(form = ~ 1|time),
#'                     data = dL, method = "ML")
#' getVarCov2(e4.gls)$Omega
#'
#' #### lvm model ####
#' m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
#' latent(m) <- ~eta
#' e <- estimate(m, dW)
#' getVarCov2(e)
#' 
#' @concept extractor
#' 
#' @export
`getVarCov2` <-
    function(object, ...) UseMethod("getVarCov2")

## * getVarCov2.lm
#' @rdname getVarCov2
#' @export
getVarCov2.lm <- function(object, ...){

    if(!inherits(object,"sCorrect")){
        .convert_sCorrect(object) <- FALSE
    }
    return(object$sCorrect$Omega)

}

## * getVarCov2.gls
#' @rdname getVarCov2
#' @export
getVarCov2.gls <- getVarCov2.lm

## * getVarCov2.lme
#' @rdname getVarCov2
#' @export
getVarCov2.lme <- getVarCov2.lm

## * getVarCov2.lvmfit
#' @rdname getVarCov2
#' @export
getVarCov2.lvmfit <- getVarCov2.lm

## * .getVarCov2 [helper]
#' @rdname getVarCov2-internal
.getVarCov2 <- function(object, param, attr.param,
                        name.endogenous, n.endogenous, ref.group, ...){

    ## ** Compute Omega
    if(n.latent>0){
        Omega <- object$conditionalMoment$value$tLambda.tiIB.Psi.iIB %*% object$conditionalMoment$value$Lambda + object$conditionalMoment$value$Sigma
    }else{
        Omega <- object$conditionalMoment$value$Sigma
    }
    
    ## ** Extract information
    var.coef <- param[attr.param$var.coef]
    cor.coef <- param[attr.param$cor.coef]

    ## ** Diagonal terms
    if(length(ref.group)>0){        
        factor.varcoef <- setNames(c(1,var.coef[-1]),
                                   attr(object$modelStruct$varStruct,"groupNames"))
        sigma2.base <- var.coef["sigma2"] * factor.varcoef[ref.group]
    }else{
        sigma2.base <- rep(var.coef["sigma2"], n.endogenous)
    }
    ## re-order according to the order of the correlation coefficients (if possible)
    if(length(cor.coef)>1 & length(var.coef)>1){
        cor.level <- gsub("corCoef","",names(cor.coef))
        var.level <- names(sigma2.base)
        var.permlevel <- .allPermutations(var.level)

        M.try <- apply(var.permlevel, MARGIN = 1, function(iLevel){
            all(apply(utils::combn(iLevel, m = 2), MARGIN = 2, FUN = paste, collapse = "") == cor.level)
        })
        if(any(M.try)){
            sigma2.base <- sigma2.base[var.permlevel[which(M.try),]]
        }
    }
    
    Omega <- diag(as.double(sigma2.base),
                  nrow = n.endogenous, ncol = n.endogenous)

    ## ** Extra-diagonal terms
    if(length(cor.coef)>0){
        index.lower <- which(lower.tri(Omega))
        index.lower.arr <- which(lower.tri(Omega),arr.ind = TRUE)
        vec.sigma.tempo <- apply(index.lower.arr,1,function(x){prod(sqrt(sigma2.base[x]))})        
        Omega[index.lower] <- cor.coef*vec.sigma.tempo
        Omega <- symmetrize(Omega, update.upper = TRUE)
    }    

    ## ** names
    if(all(!duplicated(names(sigma2.base)))){
        dimnames(Omega) <- list(names(sigma2.base), names(sigma2.base))
    }else{
        dimnames(Omega) <- list(name.endogenous, name.endogenous)
    }

    ## ** add contribution of the random effect
    if(inherits(object,"lme")){
        ran.coef <- param[attr.param$ran.coef]
        Omega <- Omega + ran.coef
    }

    
    ## ** export
    return(Omega)
}




######################################################################
### sCorrect-getVarCov2.R ends here
