### sCorrect-getVarCov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (11:00) 
## Version: 
## Last-Updated: jan 15 2020 (14:46) 
##           By: Brice Ozenne
##     Update #: 61
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
    function(fobject, param, data, ssc) UseMethod("getVarCov2")

## * getVarCov2.lm
#' @rdname getVarCov2
#' @export
getVarCov2.lm <- function(object, param = NULL, data = NULL,
                          ssc = lava.options()$ssc){

    if(is.null(object$sCorrect) || !is.null(param) || !is.null(data) || !identical(object$sCorrect$ssc$type,ssc)){
        object <- sCorrect(object, param = param, data = data, first.order = !is.null(ssc), ssc = ssc, df = NA)
    }
    Omega <- object$sCorrect$moment$Omega
    attr(Omega, "detail") <- NULL
    return(Omega)

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

## * getVarCov2.sCorrect
#' @rdname getVarCov2
getVarCov2.sCorrect <- function(object, param = NULL, data = NULL,
                               ssc = object$sCorrect$ssc$type){
    class(object) <- setdiff(class(object),"sCorrect")
    return(getVarCov2(object, param = param, data = data, ssc = ssc))
}

## * getVarCov.sCorrect
#' @rdname getVarCov2
getVarCov.sCorrect <- getVarCov2.sCorrect



######################################################################
### sCorrect-getVarCov2.R ends here